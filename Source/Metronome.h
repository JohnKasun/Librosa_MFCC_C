/*
  ==============================================================================

    Metronome.h
    Created: 9 Oct 2020 9:37:58am
    Author:  JohnK

  ==============================================================================
*/
#include <JuceHeader.h>
#pragma once
class Metronome 
{
public:

    Metronome();
    ~Metronome();
    void prepareToPlay(int samplesPerBlockExpected, double sampleRate);
    void getNextAudioBlock(const juce::AudioSourceChannelInfo& bufferToFill);
    void reset();
    void setBpm(double bpm);

    bool onMet{ false };
    bool errorMet{ false };
private:

    double mSampleRate { 0 };
    double mBpm        { 120.0 };
    int mTotalSamples  { 0 };
    int mSampleInterval{ 0 };

    juce::AudioFormatManager mFormatManager;
    std::unique_ptr <juce::AudioFormatReaderSource> pMet{ nullptr };

};