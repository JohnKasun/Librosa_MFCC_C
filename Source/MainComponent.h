#pragma once

#include <JuceHeader.h>
#include "recordAudio.h"
#include "onsetDetection.h"
#include "classifyAudio.h"
#include "Metronome.h"



class MainComponent : public juce::Component,
                      public juce::ActionListener
                      
{
public:

    MainComponent();
    ~MainComponent() override;
    void paint(juce::Graphics& g) override;
    void resized() override;
    void actionListenerCallback(const juce::String& message) override;

private:
    std::unique_ptr<juce::DialogWindow> menuAudioSetup;
    bool metEnable{ false };
    bool errorState{ false };
    bool errorMidiState{ true };
    bool doneAnalyzing{ false };

    juce::TextButton buttonRecord    { "Record" };
    juce::TextButton buttonReset     { "Reset" };
    juce::TextButton buttonStop      { "Stop Recording" };
    juce::TextButton buttonDump      { "Dump to\n .csv file" };
    juce::TextButton buttonMet       { "Metronome Off" };
    juce::TextButton buttonAnalyze   { "Analyze" };
    juce::TextButton buttonPlayMidi  { "<Set MIDI Output>" };
    juce::TextButton buttonStopMidi  { "Stop" };
    juce::TextButton audioSetupButton{ "Settings" };

    juce::TextEditor errorBox;

    juce::Slider barCount;
    juce::Label numberOfBarsLabel;
    
    juce::Slider newBpm;
    juce::Label newBpmLabel;

    //Arrays for Midi Buffer
    juce::Array<int> midiArray;
    juce::Array<int> velocityConstant;
    juce::Array<int> onsetArray;
    
    void start();
    void reset();
    void stop();
    void sendBufferVals();
    void metPressed();
    void done();
    void error();
    void dumpDataToCSV();
    void buttonAnalyzePressed();
    void buttonPlayMidiPressed();
    void ioSetup();
    
    recordAudio recorder;
    onsetDetection segments;
    classifyAudio classification;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MainComponent)
};

