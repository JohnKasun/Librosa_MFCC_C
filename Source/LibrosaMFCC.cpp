/*
  ==============================================================================

    LibrosaMFCC.cpp
    Created: 8 Dec 2020 5:08:20pm
    Author:  JohnK

  ==============================================================================
*/

#include "LibrosaMFCC.h"

Lib_Mfcc::Lib_Mfcc() : forwardFFT(fftOrder) {}
Lib_Mfcc::~Lib_Mfcc() {}

std::vector<std::vector<float>> Lib_Mfcc::doMfcc(std::vector<float> y, int sampleRate, int n_mfcc, int dct_type, int hopLength, bool centered)
{
    std::vector<std::vector<float>> error;
    if (hopLength < fftSize)
    {
        auto y_pad = y;
        if (centered)
            y_pad = padAudio(y);
        else if (y.size() < fftSize)
            return error;
        
        auto mel_basis = getMelFilterBank(sampleRate);
        auto fft = doFFT(y_pad, hopLength);
        auto signal_power = signalPower(fft);
        auto audio_filtered = doFilter(signal_power, mel_basis);
        auto cepCoeff = doDCT(audio_filtered, n_mfcc, dct_type);

        return cepCoeff;
    }
    else
        return error;
    
}

double Lib_Mfcc::freqToMel(double freq)
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

double Lib_Mfcc::melToFreq(double mel)
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

std::vector<double> Lib_Mfcc::linspace(double start_in, double end_in, int num_in)
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

std::vector<double> Lib_Mfcc::arange(double start_in, double end_in, double spacing)
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

std::vector<std::vector<float>> Lib_Mfcc::normalize(std::vector<std::vector<float>> weights, std::vector<double> mel_f)
{
    std::vector<float> enorm(melFilterNum, 0);
    auto normWeights = std::vector<std::vector<float>>(weights.size(), std::vector<float>(weights[0].size()));
    for (auto i = 0; i < melFilterNum; i++) {
        enorm[i] = (float)(2.0 / (mel_f[i + 2] - mel_f[i]));
    }

    for (auto i = 0; i < weights.size(); i++) {
        for (auto j = 0; j < weights[i].size(); j++) {
            normWeights[i][j] = weights[i][j] * enorm[i];
        }
    }

    return normWeights;
}

std::vector<std::vector<float>> Lib_Mfcc::dotProduct(std::vector<std::vector<float>> matrix1, std::vector<std::vector<float>> matrix2)
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

std::vector<float> Lib_Mfcc::padAudio(std::vector<float> y)
{
    int numPad = int(fftSize / 2);
    std::vector<float> y_pad((size_t)(y.size() + (2.0 * numPad)), 0);

    for (auto i = 0; i < y.size(); i++)
    {
        y_pad[i + numPad] = y[i];
    }

    for (auto i = 0; i < numPad; i++)
    {
        y_pad[(numPad - 1) - i] = y_pad[(numPad + 1) + i];
        y_pad[(numPad + y.size()) + i] = y_pad[(numPad + y.size() - 2) - i];
    }
    
    return y_pad;
}

std::vector<std::vector<float>> Lib_Mfcc::getMelFilterBank(double sampleRate)
{
    auto fmin_mel = freqToMel(0.0);
    auto fmax_mel = freqToMel(sampleRate / 2.0);

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

std::vector<std::vector<float>> Lib_Mfcc::doFFT(std::vector<float> audio, int hopLength)
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

std::vector<std::vector<float>> Lib_Mfcc::signalPower(std::vector<std::vector<float>> fftData)
{
    auto power = std::vector<std::vector<float>>(fftData.size(), std::vector<float>(fftData[0].size()));

    for (auto i = 0; i < fftData.size(); i++) {
        for (auto j = 0; j < fftData[0].size(); j++) {
            power[i][j] = pow(fftData[i][j], 2);
        }
    }

    return power;
}

std::vector<std::vector<float>> Lib_Mfcc::doFilter(std::vector<std::vector<float>> signal_power, std::vector<std::vector<float>> mel_basis)
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

std::vector<std::vector<float>> Lib_Mfcc::doDCT(std::vector<std::vector<float>> signal_filtered, int n_mfcc, int dct_type)
{
    std::vector<std::vector<float>> dct;
    if (dct_type == 1)
        dct = dct;
    else if (dct_type == 2)
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
        dct = result;
    }
    else
        dct = dct;

    std::vector<std::vector<float>> output(n_mfcc, std::vector<float>(dct[0].size()));
    for (auto i = 0; i < n_mfcc; i++)
    {
        for (auto j = 0; j < dct[0].size(); j++)
        {
            output[i][j] = dct[i][j];
        }
    }
    return output;
    
}