/*
  ==============================================================================

    LibrosaMFCC.h
    Created: 8 Dec 2020 5:08:20pm
    Author:  JohnK

  ==============================================================================
*/
#pragma once
#include <complex>
#include <cmath>
#include <vector>

class Lib_Mfcc
{
    static constexpr auto melFilterNum = 128;

public:
    Lib_Mfcc();
    ~Lib_Mfcc();

    //master function
    std::vector<std::vector<float>> doMfcc(std::vector<float> y, int sampleRate = 22050, int n_mfcc = 20, int dct_type = 2, bool ortho = true, int hopLength = 512, int fftSize = 2048, bool centered = true);

private:

    //helper functions
    double freqToMel(double freq);
    double melToFreq(double mel);
    std::vector<double> linspace(double start_in, double end_in, int num_in);
    std::vector<double> arange(double start_in, double end_in, double spacing);
    std::vector<std::vector<float>> normalize(std::vector<std::vector<float>> weights, std::vector<double> mel_f);
    std::vector<std::vector<float>> dotProduct(std::vector<std::vector<float>> matrix1, std::vector<std::vector<float>> matrix2);

    //process functions
    std::vector<float> padAudio(std::vector<float> y, int fftSize);
    std::vector<std::vector<float>> getMelFilterBank(double sampleRate, int fftSize);
    std::vector<std::vector<float>> doFFT(std::vector<float> audio, int hopLength, int fftSize);
    std::vector<std::complex<float>> FFT_recursion(std::vector<float> audio);
    std::vector<std::vector<float>> doFilter(std::vector<std::vector<float>> fft, std::vector<std::vector<float>> mel_basis);
    std::vector<std::vector<float>> doDCT(std::vector<std::vector<float>> signal_filtered, int n_mfcc, int dct_type, bool ortho);

    double pi = 3.14159265358;
};

