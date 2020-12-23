/*
  ==============================================================================

    recordAudio.cpp
    Created: 2 Oct 2020 12:47:40pm
    Author:  JohnK

  ==============================================================================
*/
#include <iostream>
#include "recordAudio.h"


recordAudio::recordAudio()
{
    setAudioChannels(1, 2);

    //inherited JUCE class that handles changing/initializing audio devices
    audioSetupComp.reset(new juce::AudioDeviceSelectorComponent(deviceManager, 1, 1, 2, 2, false, true, false, true));
    audioSetupComp->setSize(600, 400);

    //Checks if audio devices are configured correctly on startup
    auto* device = deviceManager.getCurrentAudioDevice();

    auto activeInputChannels = device->getActiveInputChannels();
    numInputChannels = activeInputChannels.getHighestBit() + 1;

    auto activeOutputChannels = device->getActiveOutputChannels();
    numOutputChannels = activeOutputChannels.getHighestBit() + 1;

    auto sampleRate = device->getCurrentSampleRate();

    if ((device == nullptr) || (numInputChannels == 0) || (numOutputChannels != 2))
    {
        errored = true;
    }
    else
    {

        errored = false;
        
        if (sampleRate == 44100 || sampleRate == 48000)
            erroredSR = false;
        else
            erroredSR = true;
           
    }

    if (metronome.errorMet)
        erroredMet = true;
}

recordAudio::~recordAudio()
{
    bufferRecordedAudio.setSize(0,0);
    bufferMidi.clear();

    shutdownAudio();
}



void recordAudio::getNextAudioBlock(const juce::AudioSourceChannelInfo& bufferToFill)
{

    if (isRecording)
    {

        mBufferSize = bufferToFill.numSamples;

        //Enters this block if the number of samples left to record is less than the buffer length
        if ((startSample + mBufferSize) >= bufferRecordedAudio.getNumSamples())
        {
            auto remainder = bufferRecordedAudio.getNumSamples() - startSample;
            
            auto* readerInput  = bufferToFill.buffer->getReadPointer(0);
            auto* writerRecord = bufferRecordedAudio .getWritePointer(0, startSample);

                for (auto sample = 0; sample < remainder; ++sample)
                {
                    writerRecord[sample] = readerInput[sample];
                }
            
            startSample += remainder;
            stopRecording();
        }
        //Enters this block if the entire buffer will be recorded
        else
        {
            auto* readerInput = bufferToFill.buffer->getReadPointer(0);
            auto* writerRecord = bufferRecordedAudio.getWritePointer(0, startSample);
                
                for (auto sample = 0; sample < mBufferSize; ++sample)
                {
                    writerRecord[sample] = readerInput[sample];
                }
            
            startSample += mBufferSize;
        }

        bufferToFill.buffer->clear(0, mBufferSize);
        metronome.getNextAudioBlock(bufferToFill);
        
    }
    else
    {
        bufferToFill.clearActiveBufferRegion();
    }
}

void recordAudio::prepareToPlay(int samplesPerBlockExpected, double sampleRate)
{
    auto* device = deviceManager.getCurrentAudioDevice();
 
    auto activeInputChannels = device->getActiveInputChannels();
    numInputChannels = activeInputChannels.getHighestBit() + 1;

    auto activeOutputChannels = device->getActiveOutputChannels();
    numOutputChannels = activeOutputChannels.getHighestBit() + 1;

    //Ensures selected audio device has 1 input and 2 outputs
    if ((device == nullptr) || (numInputChannels == 0) || (numOutputChannels != 2))
    {
        sendActionMessage("Device Error");
    }
    else
    {
 
        if (mSampleRate != sampleRate)
        {
            sendActionMessage("Sample Rate Change");
        }

        mSampleRate = sampleRate;
        mBufferSize = samplesPerBlockExpected;

        if (mSampleRate == 44100 || mSampleRate == 48000)
            sendActionMessage("Device Success");
        else
            sendActionMessage("Sample Rate Error");
            

        metronome.prepareToPlay(mBufferSize, mSampleRate);

        createAudioBuffer(mBar, mBpm);
    }

    //Initializes Midi Output
    midiOutput = deviceManager.getDefaultMidiOutput();
    if (midiOutput != nullptr)
    {
        midiOutput->startBackgroundThread();
        sendActionMessage("Midi Success");
    }
    else
        sendActionMessage("Midi Error");
    
}

void recordAudio::releaseResources()
{}


void recordAudio::startRecording()
{
    isRecording = true;
}

void recordAudio::stopRecording()
{
    isRecording = false;
    sendActionMessage("Done Recording");

    //Changes length of buffer after recording has stopped.  Only makes a difference if recording was stopped prematurely
    bufferRecordedAudio.setSize(numInputChannels, startSample, true, true, false);
}

void recordAudio::resetRecording()
{

    isRecording = false;
    startSample = 0;
    metronome.reset();

    createAudioBuffer(mBar, mBpm);
}

void recordAudio::createAudioBuffer(double numBar, double bpm)
{
    mBar = numBar;
    mBpm = bpm;
    metronome.setBpm(mBpm);

    auto samplesToAllocate = mBar * (60.0 / mBpm) * mSampleRate * 4;
    bufferRecordedAudio.setSize(numInputChannels, (int)samplesToAllocate);
    bufferRecordedAudio.clear();
}

void recordAudio::metEnabled(bool enable)
{
    metronome.onMet = enable;
}

void recordAudio::fillMidiBuffer(std::vector<int> onsetVec, std::vector <int> drumVec, std::vector<int> velVec)
{
    bufferMidi.clear();
    bufferMidi.ensureSize(onsetVec.size());

    //creates Midi notes for each onset and adds them into a Midi buffer
    for (auto i = 0; i < onsetVec.size(); ++i)
    {
        auto message = juce::MidiMessage::noteOn(1, drumVec[i], (juce::uint8) velVec[i]);
        bufferMidi.addEvent(message, onsetVec[i]);

        auto messageOff = juce::MidiMessage::noteOff(message.getChannel(), message.getNoteNumber());
        bufferMidi.addEvent(messageOff, onsetVec[i] + 5000);
    }
    sendActionMessage("Done Analyzing");
}

void recordAudio::doAlgorithm()
{
    auto noveltyFunction = onset.makeNoveltyFunction(bufferRecordedAudio, bufferRecordedAudio.getNumSamples(), mSampleRate);
    auto peaks = onset.pickPeaks(noveltyFunction);
    auto peaksInSeconds = onset.convertIndiciesToTime(peaks);
    auto peaksInSamples = onset.convertTimeToSamples(peaksInSeconds);

    std::vector<int> onsetVec = peaksInSamples;

    std::vector<int> drumVec = classification.splitAudio(bufferRecordedAudio, onsetVec, mSampleRate);

    std::vector<int> velVec(onsetVec.size(), 100);

    fillMidiBuffer(onsetVec, drumVec, velVec);
}

void recordAudio::outputMidi()
{
    midiOutput->clearAllPendingMessages();
    midiOutput->sendBlockOfMessages(bufferMidi, juce::Time::getMillisecondCounter(), mSampleRate);
}

void recordAudio::stopOutputMidi()
{
    midiOutput->clearAllPendingMessages();
}

void recordAudio::testAlgorithm()
{
    /*tester();

    auto noveltyFunction = onset.makeNoveltyFunction(bufferRecordedAudio, bufferRecordedAudio.getNumSamples(), mSampleRate);
    auto peaks = onset.pickPeaks(noveltyFunction);
    auto peaksInSeconds = onset.convertIndiciesToTime(peaks);
    auto peaksInSamples = onset.convertTimeToSamples(peaksInSeconds);
    onset.testSegmentation(noveltyFunction, peaks, peaksInSeconds);*/

    classification.tester(bufferRecordedAudio, mSampleRate);
}

void recordAudio::tester()
{
    juce::File myFile;

    auto parentDir = juce::File::getSpecialLocation(juce::File::userDocumentsDirectory);

    myFile = parentDir.getChildFile("Test_Audio_Recording.csv");
    myFile.deleteFile();

    juce::FileOutputStream output2(myFile);
    output2.writeString("Channel1\n");
    
    for (auto sample = 0; sample < bufferRecordedAudio.getNumSamples(); ++sample)
    {
        
        juce::String dataString1 = (juce::String) (double)bufferRecordedAudio.getSample(0, sample);
        output2.writeString(dataString1);

        output2.writeString("\n");
    }
    output2.flush();  
    myFile.~File();
    DBG("Done writing");
}

