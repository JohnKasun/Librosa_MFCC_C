/*
  ==============================================================================

    Metronome.cpp
    Created: 9 Oct 2020 9:37:58am
    Author:  JohnK

  ==============================================================================
*/

#include "Metronome.h"
Metronome::Metronome()
{
    //This block of code finds the metronome wav file
    mFormatManager.registerBasicFormats();

    juce::File myFile{ juce::File::getSpecialLocation(juce::File::SpecialLocationType::userDocumentsDirectory) };
    auto beatvoxDir = myFile.findChildFiles(juce::File::TypesOfFileToFind::findDirectories, true, "BeatVOX_main");
    auto addFilesDir = beatvoxDir[0].findChildFiles(juce::File::TypesOfFileToFind::findDirectories, true, "AdditionalFiles");

    auto metFile = addFilesDir[0].findChildFiles(juce::File::TypesOfFileToFind::findFiles, true, "Cowbell-2.wav");

    if (metFile[0].exists())
    {
        auto formatReader = mFormatManager.createReaderFor(metFile[0]);
        pMet.reset(new juce::AudioFormatReaderSource(formatReader, true));
    }
    else
    {
        errorMet = true;
    }

}

Metronome::~Metronome()
{}

void Metronome::setBpm(double bpm)
{
    mBpm = bpm;
    mSampleInterval = int((60.0 / mBpm) * mSampleRate);
}

void Metronome::prepareToPlay(int samplesPerBlockExpected, double sampleRate)
{
    mSampleRate = sampleRate;
    mSampleInterval = int((60.0 / mBpm) * mSampleRate);

    if (pMet != nullptr)
    {
        //Initializes metronome wav file for use
        pMet->prepareToPlay(samplesPerBlockExpected, sampleRate);
    }
}

void Metronome::getNextAudioBlock(const juce::AudioSourceChannelInfo& bufferToFill)
{
  
    if (onMet)
    {
        //Makes a click when record button is first pressed
        if (mTotalSamples == 0)
        {
            pMet->setNextReadPosition(0);
            pMet->getNextAudioBlock(bufferToFill);
        }

        auto bufferSize = bufferToFill.numSamples;

        auto mSamplesToGo = mTotalSamples % mSampleInterval;

        //Enters loop only if click is supposed to occur within this audio buffer block
        if (((mSamplesToGo + bufferSize) >= mSampleInterval))
        {
            pMet->setNextReadPosition(0);

            auto timeToStartClick = mSampleInterval - mSamplesToGo - 1;
            for (auto sample = 0; sample < bufferSize; sample++)
            {
                if (sample == timeToStartClick)
                {
                    pMet->getNextAudioBlock(bufferToFill);
                }
            }

        }

        mTotalSamples += bufferSize;

    }
 
}

void Metronome::reset()
{
    mTotalSamples = 0;
}

