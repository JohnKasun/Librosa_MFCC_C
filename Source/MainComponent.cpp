#include "MainComponent.h"


//==============================================================================
MainComponent::MainComponent()
{
    setOpaque(true);

    //initalizies openable window that allows user to configure audio devices and midi output
    menuAudioSetup.reset(new juce::DialogWindow("Settings", juce::Colours::black, true, true));

    //initialize bar slider and bpm slider
    addAndMakeVisible (barCount);
    barCount.setRange(1, 16,1);
    barCount.setTextValueSuffix (" Bars");
    addAndMakeVisible (newBpm);
    newBpm.setRange(20, 200, 1);
    newBpm.setTextValueSuffix (" BPM");
    newBpm.setValue(120);
    barCount.setValue(4);
    
    
    //initialize bar slider label and bpm slider label
    addAndMakeVisible (numberOfBarsLabel);
    numberOfBarsLabel.setText ("Number of Bars", juce::dontSendNotification);
    numberOfBarsLabel.attachToComponent (&barCount, true);
    
    addAndMakeVisible (newBpmLabel);
    newBpmLabel.setText ("BPM", juce::dontSendNotification);
    newBpmLabel.attachToComponent (&newBpm, true);
    

    //initialize all buttons and sliders with lambdas
    //set value change detected from sliders
    barCount.onValueChange = [this]() { sendBufferVals(); };
    newBpm.onValueChange = [this]() { sendBufferVals(); };
    
    addAndMakeVisible(buttonRecord);
    buttonRecord.setEnabled(true);
    buttonRecord.onClick = [this]() { start(); };
    
    addAndMakeVisible(buttonReset);
    buttonReset.setEnabled(false);
    buttonReset.onClick = [this]() { reset(); };
    
    addAndMakeVisible(buttonStop);
    buttonStop.setEnabled(false);
    buttonStop.onClick = [this]() { stop(); };

    addAndMakeVisible(buttonDump);
    buttonDump.setEnabled(true);
    buttonDump.onClick = [this]() { dumpDataToCSV(); };

    addAndMakeVisible(buttonMet);
    buttonMet.setEnabled(true);
    buttonMet.setClickingTogglesState(true);
    buttonMet.setToggleState(false, juce::NotificationType::dontSendNotification);
    buttonMet.onClick = [this]() { metPressed(); };

    addAndMakeVisible(audioSetupButton);
    audioSetupButton.setEnabled(true);
    audioSetupButton.onClick = [this]() 
    { 
        stop();
        ioSetup(); 
    };

    addAndMakeVisible(errorBox);
    errorBox.setReadOnly(true);
    errorBox.setColour(juce::TextEditor::backgroundColourId, juce::Colour(0x32ffffff));
    errorBox.setColour(juce::TextEditor::outlineColourId, juce::Colour(0x1c000000));
    errorBox.setColour(juce::TextEditor::shadowColourId, juce::Colour(0x16000000));

    addAndMakeVisible(buttonAnalyze);
    buttonAnalyze.setEnabled(false);
    buttonAnalyze.onClick = [this]() { buttonAnalyzePressed(); };

    addAndMakeVisible(buttonPlayMidi);
    buttonPlayMidi.setEnabled(false);
    buttonPlayMidi.onClick = [this]() { buttonPlayMidiPressed(); };

    addAndMakeVisible(buttonStopMidi);
    buttonStopMidi.setEnabled(false);
    buttonStopMidi.onClick = [this]() { recorder.stopOutputMidi(); };

    //Checks if audio devices are properly configured on startup
    if (recorder.errored)
    {
        error();
        errorBox.setText("ERROR -- RECONFIGURE AUDIO DEVICES");
    }
    else
    {
        if (recorder.erroredSR)
        {
            error();
            errorBox.setText("Please choose a different sample rate");
        }
        else
        {
            errorBox.setText("All Devices Successfully configured");
            sendBufferVals();
        }

    }

    if (recorder.erroredMet)
    {
        buttonMet.setEnabled(false);
        errorBox.setText("ERROR -- Metronome Disabled -- Couldn't find Cowbell_2.wav");
    }

    //When user changes audio devices, the changeListenerCallback function is called
    recorder.addActionListener(this);

    setSize(450, 430);


}

MainComponent::~MainComponent()
{
    recorder.removeActionListener(this);
}

void MainComponent::paint (juce::Graphics& g)
{
    g.fillAll (getLookAndFeel().findColour (juce::ResizableWindow::backgroundColourId));

    g.setColour(juce::Colours::antiquewhite);
    g.setFont(juce::Font("Times New Roman", 20.0f, juce::Font::bold));
    g.drawText("BeatVOX", getWidth()/4, 10, getWidth()/2, 30, juce::Justification::centred, true);
}

void MainComponent::resized()
{
    buttonRecord    .setBounds(10,  50, getWidth() - 20, 30);
    buttonStop      .setBounds(10,  90, getWidth() - 20, 30);
    buttonReset     .setBounds(10, 130, getWidth() - 20, 30);
    buttonDump      .setBounds(10, 170, getWidth() - 20, 30);
    barCount        .setBounds(10, 210, getWidth() - 20, 30);
    newBpm          .setBounds(10, 250, getWidth() - 20, 30);
    buttonMet       .setBounds(10, 290, getWidth() /  4, 30);
    buttonAnalyze   .setBounds(getWidth()/4, 330, getWidth()/2, 30);
    buttonPlayMidi  .setBounds(getWidth()/4, 370, getWidth()/4, 20);
    buttonStopMidi  .setBounds(getWidth()/2, 370, getWidth()/4, 20);
    errorBox        .setBounds(10, 405, getWidth() - 20, 20);
    audioSetupButton.setBounds(360, 15, getWidth() /  8, 20);
}

void MainComponent::start()
{
   
    recorder.startRecording();

    buttonRecord.setButtonText("Recording...");
    buttonRecord.setColour(juce::TextButton::buttonColourId, juce::Colours::red);

    buttonRecord.setEnabled(false);
    buttonReset .setEnabled(false);
    buttonStop  .setEnabled(true);
    barCount    .setEnabled(false);
    newBpm      .setEnabled(false);
    buttonMet   .setEnabled(false);

}

void MainComponent::reset()
{
    recorder.resetRecording();

    buttonRecord.setButtonText("Record");
    buttonRecord.setColour(juce::TextButton::buttonColourId, getLookAndFeel().findColour(juce::TextButton::buttonColourId));
    buttonAnalyze.setButtonText("Analyze");
    buttonAnalyze.setColour(juce::TextButton::buttonColourId, getLookAndFeel().findColour(juce::TextButton::buttonColourId));

    if (errorState)
        buttonRecord.setEnabled(false);
    else
        buttonRecord.setEnabled(true);

    buttonReset .setEnabled(false);
    buttonDump  .setEnabled(false);
    barCount    .setEnabled(true);
    newBpm      .setEnabled(true);

    if (!recorder.erroredMet)
        buttonMet   .setEnabled(true);

    buttonAnalyze.setEnabled(false);
    buttonPlayMidi.setEnabled(false);
    buttonStopMidi.setEnabled(false);

    doneAnalyzing = false;

}

void MainComponent::stop()
{
    if (recorder.isRecording)
        recorder.stopRecording();
}

void MainComponent::dumpDataToCSV()
{
    recorder.testAlgorithm();
}

//graphically informs the user that recording has finished
void MainComponent::done()
{

    buttonRecord.setColour(juce::TextButton::buttonColourId, juce::Colours::green);
    buttonRecord.setButtonText("Done!");

    buttonDump .setEnabled(true);
    buttonReset.setEnabled(true);
    buttonStop .setEnabled(false);

    if (!doneAnalyzing)
        buttonAnalyze.setEnabled(true);
    
}

void MainComponent::sendBufferVals()
{
    recorder.createAudioBuffer(barCount.getValue(), newBpm.getValue());
}

void MainComponent::metPressed()
{
    if (metEnable)
    {
        metEnable = false;
        buttonMet.setButtonText("Metronome Off");
    }
    else
    {
        metEnable = true;
        buttonMet.setButtonText("Metronome On");
    }
    recorder.metEnabled(metEnable);
}

void MainComponent::buttonAnalyzePressed()
{
    recorder.doAlgorithm();
    buttonAnalyze.setEnabled(false);
}

void MainComponent::buttonPlayMidiPressed()
{
    recorder.outputMidi();
}

void MainComponent::ioSetup()
{
    menuAudioSetup->showDialog("Settings", recorder.audioSetupComp.get(), this, juce::Colours::black, true, false, false);
}

//If there is an audio device error, the user cannot perform any functions until they reconfigure their devices in settings
void MainComponent::error()
{
    buttonRecord  .setEnabled(false);
    buttonMet     .setEnabled(false);
    barCount      .setEnabled(false);
    newBpm        .setEnabled(false);
}

//This function is called when recording finishes and when algorithm analysis finishes
void MainComponent::actionListenerCallback(const juce::String& message)
{
    if (message == "Device Success")
    {
        if (errorState)
        {
            errorState = false;
            reset();
        }
 
        sendBufferVals();
        errorBox.setText("All Devices Successfully configured");
    }
    else if (message == "Midi Success")
    {
        errorMidiState = false;
        buttonPlayMidi.setButtonText("Play");

        if (doneAnalyzing)
        {
            buttonPlayMidi.setEnabled(true);
            buttonStopMidi.setEnabled(true);
        }
           
    }
    else if (message == "Device Error")
    {
        error();
        errorState = true;
        errorBox.setText("ERROR -- RECONFIGURE AUDIO DEVICES");
    }
    else if (message == "Midi Error")
    {
        errorMidiState = true;

        buttonPlayMidi.setButtonText("<Set Midi Output>");
        buttonPlayMidi.setEnabled(false);
        buttonStopMidi.setEnabled(false);

    }
    else if (message == "Sample Rate Change")
    {
        reset();
    }
    else if (message == "Sample Rate Error")
    {
        error();
        buttonAnalyze.setEnabled(false);
        buttonPlayMidi.setEnabled(false);
        buttonStopMidi.setEnabled(false);
        errorState = true;
        errorBox.setText("Please choose a different sample rate");
    }
    else if (message == "Done Recording")
    {
        done();
    }
    else if (message == "Done Analyzing")
    {
        doneAnalyzing = true;

        buttonAnalyze.setColour(juce::TextButton::buttonColourId, juce::Colours::green);
        buttonAnalyze.setButtonText("Done!");

        if (!errorMidiState)
        {
            buttonPlayMidi.setEnabled(true);
            buttonStopMidi.setEnabled(true);
        }

    }
}
     
    


