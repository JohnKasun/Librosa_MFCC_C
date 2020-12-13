/*
  ==============================================================================

    classifyAudio.cpp
    Created: 11 Nov 2020 8:04:20pm
    Author:  JohnK

  ==============================================================================
*/

#include "classifyAudio.h"
#include <math.h>
#include "svm.h"
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))


classifyAudio::classifyAudio() : forwardFFT(fftOrder)
{
    mFormatManager.registerBasicFormats();

    juce::File myFile{ juce::File::getSpecialLocation(juce::File::SpecialLocationType::userDocumentsDirectory) };
    auto mySamples = myFile.findChildFiles(juce::File::TypesOfFileToFind::findDirectories, true, "BeatVOX_main");
    auto addFiles = mySamples[0].findChildFiles(juce::File::TypesOfFileToFind::findDirectories, true, "AdditionalFiles");

    auto model1Dir = addFiles[0].findChildFiles(juce::File::TypesOfFileToFind::findFiles, true, "model_file44100");
    auto model1 = model1Dir[0].getFullPathName();
    auto model1STD = model1.toStdString();
    const char* model1Char = model1STD.c_str();

    auto model2Dir = addFiles[0].findChildFiles(juce::File::TypesOfFileToFind::findFiles, true, "model_file48000");
    auto model2 = model2Dir[0].getFullPathName();
    auto model2STD = model2.toStdString();
    const char* model2Char = model2STD.c_str();

    model44100 = svm_load_model(model1Char);
    model48000 = svm_load_model(model2Char);

};

classifyAudio::~classifyAudio() {};

void classifyAudio::tester(juce::AudioBuffer<float> buffer, double sampleRate)
{

    juce::File myFile{ juce::File::getSpecialLocation(juce::File::SpecialLocationType::userDocumentsDirectory) };
    auto mySamples = myFile.findChildFiles(juce::File::TypesOfFileToFind::findFiles, true, "impulse-response.wav");

    auto reader = mFormatManager.createReaderFor(mySamples[0]);
    juce::AudioSampleBuffer bufferTest((int)reader->numChannels, (int)reader->lengthInSamples);
    reader->read(&bufferTest, 0, (int)reader->lengthInSamples, 0, true, true);

    std::vector<float>audio(bufferTest.getNumSamples(), 0);
    for (int i = 0; i < bufferTest.getNumSamples(); i++) {
        audio[i] = bufferTest.getSample(0, i);
    }
   
    mSampleRate = 48000;

    float absmax = 0;
    for (size_t i = 0; i < audio.size(); ++i)
    {
        absmax = std::max(absmax, std::abs(audio[i]));
    }
    for (size_t i = 0; i < audio.size(); ++i)
    {
        audio[i] = audio[i] / absmax;
    }

    auto mfcc = Mfcc.doMfcc(audio, 48000, 20, 2, true, 512, true);

    //testAccuracy1D(audio);
    testAccuracy2D(mfcc);

}

std::vector<int> classifyAudio::splitAudio(juce::AudioBuffer<float>buffer, std::vector<int>peaks, double sampleRate)
{
    mSampleRate = sampleRate;

    juce::File myFile{ juce::File::getSpecialLocation(juce::File::SpecialLocationType::userDocumentsDirectory) };
    auto mySamples = myFile.findChildFiles(juce::File::TypesOfFileToFind::findFiles, true, "battleclip_daq.wav");

    auto reader = mFormatManager.createReaderFor(mySamples[0]);
    juce::AudioSampleBuffer bufferTest((int)reader->numChannels, (int)reader->lengthInSamples);
    reader->read(&bufferTest, 0, (int)reader->lengthInSamples, 0, true, true);

    std::vector<float>audio(bufferTest.getNumSamples(), 0);
    for (int i = 0; i < bufferTest.getNumSamples(); i++) {
       audio[i] = bufferTest.getSample(0, i);
    }

    testAccuracy1D(audio);

    /*std::vector<float>audio(buffer.getNumSamples(), 0);
    for (int i = 0; i < buffer.getNumSamples(); i++) {
        audio[i] = buffer.getSample(0, i);
    }*/

    float absmax = 0;
    for (size_t i = 0; i < audio.size(); ++i)
    {
        absmax = std::max(absmax, std::abs(audio[i]));
    }
    for (size_t i = 0; i < audio.size(); ++i)
    {
        audio[i] = audio[i] / absmax;
    }
    

    auto mel_basis = getMelFilterBank(mSampleRate);
  
    //initialize midi drum array
    std::vector<int> drumArray(peaks.size(),0);
  
    for (auto i = 0; i < peaks.size(); i++)
    {
        int length = int(mSampleRate * 0.04);

        auto start_ind = peaks[i];
        auto end_ind = start_ind + length;

        if (end_ind > audio.size())
        {
            end_ind = (int)audio.size();
        }

        if (i != (peaks.size() - 1))
        {
            if (end_ind >= peaks[i+1])
                end_ind = peaks[i+1];
        }

        length = end_ind - start_ind;

        int numPad = int(fftSize / 2);

        std::vector<float> sectionTest(length, 0);
        for (auto j = 0; j < length; j++)
            sectionTest[j] = audio[start_ind + j];


        std::vector<float> section((size_t)(length + (2.0 * numPad)), 0);
   
        for (auto j = 0; j < length; j++)
        {
            section[j + numPad] = audio[j + start_ind];
        }
     
        for (auto j = 0; j < numPad; j++)
        {
            section[(numPad - 1) - j] = section[(numPad + 1) + j];
            section[(numPad + length) + j] = section[(numPad + length - 2) - j];
        }

        auto fft = doFFT(section);
        auto signal_power = signalPower(fft);
        auto audio_filtered = doFilter(signal_power, mel_basis);
        auto cepCoeff = constructDCT(audio_filtered);
        auto meanCepCoeff = meanMfcc(cepCoeff);
        auto normalizedVec = normalizeFeatures(meanCepCoeff);
 
        struct svm_node* testnode = Malloc(svm_node, normalizedVec.size() + 1);

        for (int j = 0; j < normalizedVec.size(); j++)
        {
            testnode[j].index = j+1;
            testnode[j].value = normalizedVec[j];
        }
        testnode[normalizedVec.size()].index = -1;

        auto label = 0.0;
        if (mSampleRate == 48000)
        {
            label = svm_predict(model48000, testnode);
            drumArray[i] = label;
        }
        else
        {
             label = svm_predict(model44100, testnode);
             drumArray[i] = label;
        }
        
    }
    return drumArray;
}

std::vector<std::vector<float>> classifyAudio::doFFT(std::vector<float> audio)
{
    int numOfFFTs = 1 + int((audio.size() - fftSize) / hopLength);
    std::vector<std::vector<float>> fftData(numOfFFTs, std::vector<float>(fftSize / 2 + 1));

    for (int i = 0; i < numOfFFTs; i++) {
    std::vector<float> audioData(fftSize * 2);

        // Prepare for FFT
        for (int n = 0; n < (fftSize); n++) {
            float hannWindowMultiplier = (float)(0.5 * (1.0 - cos(2.0 * pi * n / ((float)fftSize))));
            audioData[n] = hannWindowMultiplier * (float)audio[n + (hopLength * i)];
        }

        // JUCE FFT
        forwardFFT.performFrequencyOnlyForwardTransform(audioData.data());

        // Take only positive frequency part of fft
        std::vector<float>posfftData(1 + (fftSize / 2), 0);
        
        for (int j = 0; j <= fftSize / 2; j++) 
        {
             posfftData[j] = audioData[j];
        }
       
        fftData[i] = posfftData;
    }

    return fftData;
}

std::vector<std::vector<float>> classifyAudio::normalize(std::vector<std::vector<float>> weights, std::vector<double> mel_f)
{
    std::vector<float> enorm(melFilterNum, 0);
    auto normWeights = std::vector<std::vector<float>>(weights.size(), std::vector<float>(weights[0].size()));
    for (auto i = 0; i < melFilterNum; i++) {
        enorm[i] = (float)(2.0 / (mel_f[i + 2] - mel_f[i]));
    }

    for (auto i = 0; i < weights.size(); i++) {
        for(auto j = 0; j < weights[i].size(); j++) {
            normWeights[i][j] = weights[i][j] * enorm[i];
        }
    }

    return normWeights;
}

std::vector<std::vector<float>> classifyAudio::signalPower(std::vector<std::vector<float>> fftData)
{
    auto power = std::vector<std::vector<float>>(fftData.size(), std::vector<float>(fftData[0].size()));

    for (auto i = 0; i < fftData.size() ; i++) {
        for (auto j = 0; j < fftData[0].size() ; j++) {
            power[i][j] = pow(fftData[i][j], 2);
        }
    }

    return power;
}

double classifyAudio::freqToMel(double freq)
{

    auto f_min = 0.0;
    auto f_sp = 200.0 / 3;

    auto mels = (freq - f_min) / f_sp;

    auto min_log_hz = 1000.0;
    auto min_log_mel = (min_log_hz - f_min) / f_sp;
    auto logstep = log(6.4) / 27.0;

    if (freq >= min_log_hz)
        mels = min_log_mel + log(freq / min_log_hz) / logstep;

    return mels;
}

double classifyAudio::melToFreq(double mel)
{

    auto f_min = 0.0;
    auto f_sp = 200.0 / 3;

    auto freqs = (f_min + f_sp) * mel;

    auto min_log_hz = 1000.0;
    auto min_log_mel = (min_log_hz - f_min) / f_sp;
    auto logstep = log(6.4) / 27.0;

    if (mel >= min_log_mel)
        freqs = min_log_hz * exp(logstep * (mel - min_log_mel));

    return freqs;
}

std::vector<std::vector<float>> classifyAudio::getMelFilterBank(double sampleRate)
{
    auto fmin_mel = freqToMel(0);
    auto fmax_mel = freqToMel(sampleRate/2);

    std::vector<std::vector<float>> weights(melFilterNum, std::vector<float>(1 + fftSize / 2));

    auto fftFreqs = linspace(0, sampleRate / 2, int(1 + fftSize / 2));
    auto mel_f = linspace(fmin_mel, fmax_mel, melFilterNum + 2);

    for (auto i = 0; i < mel_f.size(); i++)
        mel_f[i] = melToFreq(mel_f[i]);

    std::vector<float> fdiff(mel_f.size() - 1, 0);

    for (auto i = 0; i < fdiff.size(); i++)
        fdiff[i] = (float)(mel_f[i + 1] - mel_f[i]);

    std::vector<std::vector<float>> ramps(mel_f.size(), std::vector<float>(fftFreqs.size()));
    for (auto i = 0; i < mel_f.size(); i++)
        for (auto j = 0; j < fftFreqs.size(); j++)
            ramps[i][j] = (float)mel_f[i] - (float)fftFreqs[j];

    for (auto i = 0; i < melFilterNum; i++)
    {
        std::vector<float>lower(ramps[0].size(), 0);
        std::vector<float>upper(ramps[0].size(), 0);
        for (auto j = 0; j < ramps[0].size(); j++)
        {
            lower[j] = -ramps[i][j] / fdiff[i];
            upper[j] = ramps[i + 2][j] / fdiff[i + 1];
        }

        std::vector<float> minimum(lower.size(), 0);
        for (auto x = 0; x < lower.size(); x++)
        {
            if (lower[x] <= upper[x])
                minimum[x] = lower[x];
            else
                minimum[x] = upper[x];
        }
        
        std::vector<float> maximum(minimum.size(), 0);
        for (auto x = 0; x < minimum.size(); x++)
        {
            if (minimum[x] > maximum[x])
                maximum[x] = minimum[x];
        }

        weights[i] = maximum;

    }

    weights = normalize(weights, mel_f);

    return weights;

}

std::vector<std::vector<float>> classifyAudio::doFilter(std::vector<std::vector<float>> signal_power, std::vector<std::vector<float>> mel_basis)
{
    std::vector<std::vector<float> > trans_vec(signal_power[0].size(), std::vector<float>(signal_power.size()));

    for (int i = 0; i < signal_power.size(); i++)
    {
        for (int j = 0; j < signal_power[i].size(); j++)
        {
            trans_vec[j][i] = signal_power[i][j];
        }
    }

    
    auto dot = dotProduct(mel_basis, trans_vec);

    auto amin = (float)1e-10;
    auto topDB = (float)80.0;

    std::vector<std::vector<float> > audio_log(dot.size(), std::vector<float>(dot[0].size()));

    for (auto i = 0; i < audio_log.size(); i++) {
        for (auto j = 0; j < audio_log[0].size(); j++) 
        {
            auto val1 = dot[i][j];
            if (amin >= dot[i][j])
                val1 = amin;

            audio_log[i][j] = (float)(10.0 * log10(val1));
        }
    }

    auto max = audio_log[0][0];
    for (auto i = 0; i < audio_log.size(); i++)
    {
        for (auto j = 0; j < audio_log[0].size(); j++)
        {
            if (audio_log[i][j] > max)
                max = audio_log[i][j];
        }
    }

    for (auto i = 0; i < audio_log.size(); i++) {
        for (auto j = 0; j < audio_log[0].size(); j++)
        {
            auto newVal = max - topDB;
            if (newVal > audio_log[i][j])
                audio_log[i][j] = newVal;
        }
    }

    return (audio_log);

}

std::vector<std::vector<float>> classifyAudio::constructDCT(std::vector<std::vector<float>> signal_filtered)
{

    auto col = signal_filtered[0].size();
    auto N = signal_filtered.size();
    std::vector<std::vector<float>> result(N, std::vector<float>(col));

    for (auto c = 0; c < col; c++)
    {
        for (auto k = 0; k < N; k++)
        {
            auto sum = 0.0;
            for (auto n = 0; n < N; n++)
            {
                if (k == 0)
                    sum += sqrt(1.0 / (4.0 * N)) * 2.0 * signal_filtered[n][c] * cos((pi * k) * (2.0 * n + 1.0) / (2.0 * N));
                else
                    sum += sqrt(1.0 / (2.0 * N)) * 2.0 * signal_filtered[n][c] * cos((pi * k) * (2.0 * n + 1.0) / (2.0 * N));

                
            }
            result[k][c] = (float)sum;
        }
    }

    std::vector<std::vector<float>> output(dctFilterNum,std::vector<float>(result[0].size()));

    for (auto i = 0; i < dctFilterNum; i++)
    {
        for (auto j = 0; j < result[0].size(); j++)
        {
            output[i][j] = result[i][j];
        }
    }
    return output;
}

std::vector<std::vector<float>> classifyAudio::dotProduct(std::vector<std::vector<float>> matrix1, std::vector<std::vector<float>> matrix2)
{
     std::vector<std::vector<float>> output(matrix1.size(), std::vector<float>(matrix2[0].size()));

    if (matrix1[0].size() != matrix2.size())
    {
        DBG("cannot perform dot product");
    }
    else
    {
        for (auto i = 0; i < matrix1.size(); i++)
        {
            for (auto j = 0; j < matrix2[0].size(); j++)
            {
                for (auto k = 0; k < matrix1[0].size(); k++)
                {
                    output[i][j] += matrix1[i][k] * matrix2[k][j];
                }
            }
        }
    }
    return output;
}

std::vector<double> classifyAudio::meanMfcc(std::vector<std::vector<float>> matrix)
{
    std::vector<double>mean(matrix.size());
    double rowsum;
    int index;

    for (int i = 0; i < matrix.size(); i++)
    {
        rowsum = 0;
        index = 0;
        for (int j = 0; j < matrix[0].size(); j++)
        {
            index += 1;
            rowsum += matrix[i][j];
        }
        mean[i] = rowsum/index;        
    }
    return mean;
}

std::vector<double> classifyAudio::normalizeFeatures(std::vector<double> featureVec)
{
    if (featureVec.size() != minVals44100.size() || featureVec.size() != maxVals44100.size())
        DBG("cannot perform feature normalization");

    std::vector<double> minVals, maxVals;

    if (mSampleRate == 48000)
    {
        minVals = minVals48000;
        maxVals = maxVals48000;
    }
    else
    {
        minVals = minVals44100;
        maxVals = maxVals44100;
    }
    std::vector<double> normFeatureVec(featureVec.size());

    for (auto i = 0; i < minVals.size(); i++)
    {
        normFeatureVec[i] = (featureVec[i] - minVals[i]) / (maxVals[i] - minVals[i]);
    }
    return normFeatureVec;
}

std::vector<double> classifyAudio::linspace(double start_in, double end_in, int num_in)
{

    std::vector<double> linspaced;

    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double num = static_cast<double>(num_in);

    if (num == 0) { return linspaced; }
    if (num == 1)
    {
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for (int i = 0; i < num - 1; ++i)
    {
        linspaced.push_back(start + delta * i);
    }
    linspaced.push_back(end); 
    return linspaced;
}

std::vector<double> classifyAudio::arange(double start_in, double end_in, double spacing)
{
    std::vector<double>aranged;
    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double num = static_cast<double>(spacing);
    double current = start;
    aranged.push_back(start);


    while ((current + num) < end)
    {
        aranged.push_back(current + spacing);
        current += spacing;
    }
    return aranged;
}

void classifyAudio::testAccuracy2D(std::vector<std::vector<float>> cepCoeff) {
    juce::File myFile;

    auto parentDir = juce::File::getSpecialLocation(juce::File::userDocumentsDirectory);

    myFile = parentDir.getChildFile("test_Mfcc.csv");
    myFile.deleteFile();

    juce::FileOutputStream output2(myFile);

    for (auto frame = 0; frame < cepCoeff[0].size(); ++frame)
    {
        output2.writeString("frame" + (juce::String)frame + ",");
    }

    output2.writeString("\n");

    for (auto mfcc = 0; mfcc < cepCoeff.size(); ++mfcc)
    {

        for (auto frame = 0; frame < cepCoeff[0].size(); ++frame)
        {

            juce::String dataString1 = (juce::String) cepCoeff[mfcc][frame];
            output2.writeString(dataString1);
            output2.writeString(",");

        }

        output2.writeString("\n");
    }
    output2.flush();
    myFile.~File();

}

void classifyAudio::testAccuracy1D(std::vector<float> section)
{
    juce::File myFile;

    auto parentDir = juce::File::getSpecialLocation(juce::File::userDocumentsDirectory);

    myFile = parentDir.getChildFile("test_Audio.csv");
    myFile.deleteFile();

    juce::FileOutputStream output2(myFile);

    output2.writeString("Data\n");

    for (auto i = 0; i < section.size(); ++i)
    {
        juce::String dataString1 = (juce::String) section[i];
        output2.writeString(dataString1);
        output2.writeString("\n");
    }
    output2.flush();
    myFile.~File();
}


