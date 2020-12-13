/*
  ==============================================================================

    recordAudio.h
    Created: 2 Oct 2020 12:47:40pm
    Author:  JohnK

  ==============================================================================
*/

#include <JuceHeader.h>
#include "Metronome.h"
#include "classifyAudio.h"
#include "onsetDetection.h"
#pragma once

class recordAudio : public juce::AudioAppComponent, 
                    public juce::ActionBroadcaster
{
public:
    recordAudio();
    ~recordAudio();
    void getNextAudioBlock(const juce::AudioSourceChannelInfo& bufferToFill) override;
    void prepareToPlay(int samplesPerBlockExpected, double sampleRate) override;
    void releaseResources();

    void startRecording();
    void stopRecording();
    void createAudioBuffer(double numBars, double bpm);
    void resetRecording();
    void metEnabled(bool enable);
    void tester();
    void testAlgorithm();
    void doAlgorithm();
    void fillMidiBuffer(std::vector<int> onsetVec, std::vector <int> drumVec, std::vector<int> velVec);
    void outputMidi();
    void stopOutputMidi();

    bool isRecording{ false };
    bool errored{ false };
    bool erroredMet{ false };
    bool erroredSR{ false };
 
    std::unique_ptr<juce::AudioDeviceSelectorComponent> audioSetupComp;
    juce::AudioBuffer<float> bufferRecordedAudio;
    double mSampleRate{ 0.0 };

private:
    
    double mBpm          { 120.0 };
    double mBar          { 4.0 };
    int startSample      { 0 };
    int numInputChannels { 1 };
    int numOutputChannels{ 2 };
    int mBufferSize      { 0 };

    juce::MidiBuffer bufferMidi;
    juce::MidiOutput* midiOutput{ nullptr };

    Metronome metronome;
    classifyAudio classification;
    onsetDetection onset;
    
};
