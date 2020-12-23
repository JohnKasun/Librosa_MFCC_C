/*
  ==============================================================================

    onsetDetection.cpp
    Created: 14 Oct 2020 2:38:45pm
    Author:  Walter Kopacz

  ==============================================================================
*/

#include <math.h>
#include "onsetDetection.h"

onsetDetection::onsetDetection() : forwardFFT(fftOrder),
hannWindow(fftSize, juce::dsp::WindowingFunction<float>::hann){

     
}

onsetDetection:: ~onsetDetection() {
    //peaks.clear();
    //peaksInSeconds.clear();
    //
}

std::vector<float> onsetDetection::makeNoveltyFunction(juce::AudioBuffer<float>buffer, int audioNumOfSamples, double sampleRate) {
    
    mSampleRate = sampleRate;
    double pi = 3.141592653589793;
    
    //zero pad
    //int numberOfZeros = fftSize - (audioNumOfSamples%(hopLength));
    // turn buffer into vector
    //std::vector<float>audio(audioNumOfSamples+numberOfZeros, 0);
    std::vector<double>audio(audioNumOfSamples, 0);
    for(int i = 0; i < audioNumOfSamples; i++){
        audio[i] = buffer.getSample(0, i);
        audio[i] = round(audio[i]*100000000)/100000000;
    }
    // Normalize the audio
    double absmax = 0;
    for (size_t i = 0; i < audioNumOfSamples; ++i)
    {
        absmax = std::max(absmax, std::abs(audio[i]));
    }
    absmax = round(absmax*100000)/100000;
    for (size_t i = 0; i < audioNumOfSamples; ++i)
    {
        audio[i] = audio[i] / absmax;
    }
    // Allocate array for FFT data
    
    //classification.testAccuracy1(audio);

    int numOfFFTs = 1 + int((audio.size() - fftSize) / hopLength);
    std::vector<std::vector<float>> fftData(numOfFFTs, std::vector<float>(fftSize / 2 + 1));
    /*int numOfFFTs = ceil((audioNumOfSamples)/(hopLength));
    fftData.resize(numOfFFTs);*/
    
    for(int i = 0; i < numOfFFTs; i++) {
        std::vector<double> audioData(fftSize*2);
        
        // Prepare for FFT
        for (int n = 0; n < (fftSize); n++){
            double hannWindowMultiplier =(double)(0.5 * (1.0 - cos(2.0 * pi * n / ((double)fftSize))));
            audioData[n] = hannWindowMultiplier * audio[n+(hopLength*i)];
        }
        // JUCE FFT
        std::vector<float> audioDataFloat(audioData.begin(), audioData.end());
        forwardFFT.performFrequencyOnlyForwardTransform(audioDataFloat.data());

        // Take only positive frequency part of fft
        std::vector<float>posfftData(1 + (fftSize / 2), 0);

        for (int j = 0; j <= fftSize / 2; j++)
        {
            posfftData[j] = audioDataFloat[j];
        }

        fftData[i] = posfftData;
    }
    
    //classification.testAccuracy(fftData);

    //Convert to Magnitude Spectrum to dB
    for (auto i = 0; i < fftData.size(); i++)
        for (auto j = 0; j < fftData[0].size(); j++)
            fftData[i][j] = 20.0 * log10f(fftData[i][j]+1e-12);

    //Duplicate the first frame so first flux will be zero
    fftData.insert(fftData.begin(), fftData[0]);
    
    // Diff each spectrum
    std::vector<std::vector<float>> flux(fftData.size()-1);
    for(int i = 0; i < flux.size(); i++){ flux[i].resize(fftSize/2 + 1); }
    
    for(int i = 0; i < fftData.size() - 1; i++){
        for(int j = 0; j < (fftData[i].size()); j++){
            flux[i][j] = fftData[i+1][j] - fftData[i][j];
            
        }
    }
    //Halfwave rectification
    for(int i = 0; i < flux.size(); i++){
        for(int j = 0; j < flux[i].size(); j++){
            if (flux[i][j] < 0) {
                flux[i][j] = 0;
            }
        }
    }
    // square for RMS
    for (int i = 0;i < flux.size(); i++) {
        for (int j = 0; j < flux[i].size(); j++) {
            flux[i][j] = pow(flux[i][j], 2);
        }
    }
  
    
   
    
    // RMS
    std::vector<float>mean(numOfFFTs);
    
    for(int i = 0; i < flux.size(); i++){
        float sum = 0.0;
        for(int j = 0; j < flux[i].size(); j++){
            sum = sum + flux[i][j];
        }
        mean[i] = sum / flux[i].size();
    }
    
    std::vector<float>RMS(mean.size());
    RMS = mean;
    for(int i = 0; i < mean.size(); i++){
        RMS[i] = sqrt(mean[i]);
    }
    // Normalize novelty function
    float absmax_nov = 0;
    for (size_t i = 0; i < RMS.size(); ++i)
    {
        absmax_nov = std::max(absmax_nov, std::abs(RMS[i]));
    }
    for (size_t i = 0; i < RMS.size(); ++i)
    {
        RMS[i] = RMS[i] / absmax_nov;
    }
    auto noveltyFunction = RMS;

    return noveltyFunction;
    
}

std::vector<int> onsetDetection::pickPeaks(std::vector<float>noveltyFunction) {

    // Get adaptive threshold using median filter ( 5 points )
    std::vector<float>adaptiveThreshold(noveltyFunction.size(), 0);
    std::vector<float>noveltyFunctionZeroPadded(noveltyFunction.size()+4, 0); //for now assuming it adds 2 zero at the beginning and end
    for(int i = 2; i < noveltyFunctionZeroPadded.size(); i++){
        if (i < noveltyFunction.size()+2){
            noveltyFunctionZeroPadded[i] = noveltyFunction[i-2];
        }
    }
    std::array<float, 5>currentBlock;
    for(int i = 2; i < noveltyFunctionZeroPadded.size()-2; i++){
        currentBlock[0] = noveltyFunctionZeroPadded[i-2];
        currentBlock[1] = noveltyFunctionZeroPadded[i-1];
        currentBlock[2] = noveltyFunctionZeroPadded[i];
        currentBlock[3] = noveltyFunctionZeroPadded[i+1];
        currentBlock[4] = noveltyFunctionZeroPadded[i+2];
        std::sort(currentBlock.begin(), currentBlock.end());
        adaptiveThreshold[i-2] = currentBlock[2] + .2;
    }
    //pick the peaks
    std::vector<int> peaks;
    for(int i = 1; i < noveltyFunction.size()-1; i++){
        if (noveltyFunction[i-1] < noveltyFunction[i] &&
            noveltyFunction[i+1] < noveltyFunction[i] &&
            noveltyFunction[i] > adaptiveThreshold[i]){
            peaks.push_back(i);
            i = i+3;
            
        }
    }

    return peaks;
}

std::vector<float> onsetDetection::convertIndiciesToTime(std::vector<int>peaksInIndicies){
    std::vector<float> peaksInSeconds;
    if (peaksInIndicies.size() > 0){
        for (int i = 0; i < peaksInIndicies.size(); i++){
            peaksInSeconds.push_back(peaksInIndicies[i] * hopLength / mSampleRate);
        }
    }
    return peaksInSeconds;
}

std::vector<int> onsetDetection::convertTimeToSamples(std::vector<float>peaksInTime){
    std::vector<int> peaksInSamples;
    if (peaksInTime.size() > 0){
        for (int i = 0; i < peaksInTime.size(); i++){
            peaksInSamples.push_back((int)(peaksInTime[i] * mSampleRate));
        }
    }
    return peaksInSamples;
}

void onsetDetection::testSegmentation(std::vector<float>noveltyFunction, std::vector<int> peaks, std::vector<float> peaksInSeconds){
    juce::File myFile;

    auto parentDir = juce::File::getSpecialLocation(juce::File::userDocumentsDirectory);

    myFile = parentDir.getChildFile("Test_Segmentation.csv");
    myFile.deleteFile();

    juce::FileOutputStream output2(myFile);
    output2.writeString("Novelty Function\n");
    //Output npvelty fucntion to CSV
    for (auto sample = 0; sample < noveltyFunction.size(); ++sample)
    {
        juce::String dataString1 = (juce::String) noveltyFunction[sample];
        output2.writeString(dataString1);
        if (sample < peaks.size()){
            juce::String dataString2 = (juce::String) peaks[sample];
            juce::String dataString3 = (juce::String) peaksInSeconds[sample];
            output2.writeString(",");
            output2.writeString(dataString2);
            output2.writeString(",");
            output2.writeString(dataString3);
        }

        output2.writeString("\n");
    }
    output2.flush();
    myFile.~File();
    DBG("Done writing");
}


    

