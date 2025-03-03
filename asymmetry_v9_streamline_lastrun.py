#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v2024.2.4),
    on Sun Mar  2 20:17:43 2025
If you publish work using this script the most relevant publication is:

    Peirce J, Gray JR, Simpson S, MacAskill M, Höchenberger R, Sogo H, Kastman E, Lindeløv JK. (2019) 
        PsychoPy2: Experiments in behavior made easy Behav Res 51: 195. 
        https://doi.org/10.3758/s13428-018-01193-y

"""

# --- Import packages ---
from psychopy import locale_setup
from psychopy import prefs
from psychopy import plugins
plugins.activatePlugins()
from psychopy import sound, gui, visual, core, data, event, logging, clock, colors, layout, hardware
from psychopy.tools import environmenttools
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER, priority)

import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import (sin, cos, tan, log, log10, pi, average,
                   sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle, choice as randchoice
import os  # handy system and path functions
import sys  # to get file system encoding

from psychopy.hardware import keyboard

# Run 'Before Experiment' code from block_start_code
import random, numpy as np
# Run 'Before Experiment' code from rand_divider_slider
from psychopy.hardware import keyboard
kb = keyboard.Keyboard()
# --- Setup global variables (available in all functions) ---
# create a device manager to handle hardware (keyboards, mice, mirophones, speakers, etc.)
deviceManager = hardware.DeviceManager()
# ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
# store info about the experiment session
psychopyVersion = '2024.2.4'
expName = 'asymmetry_v8'  # from the Builder filename that created this script
# information about this experiment
expInfo = {
    'subj': '',
    'sess_type': 'A or B or C or D',
    'date|hid': data.getDateStr(),
    'expName|hid': expName,
    'psychopyVersion|hid': psychopyVersion,
}

# --- Define some variables which will change depending on pilot mode ---
'''
To run in pilot mode, either use the run/pilot toggle in Builder, Coder and Runner, 
or run the experiment with `--pilot` as an argument. To change what pilot 
#mode does, check out the 'Pilot mode' tab in preferences.
'''
# work out from system args whether we are running in pilot mode
PILOTING = core.setPilotModeFromArgs()
# start off with values from experiment settings
_fullScr = True
_winSize = [1920, 1080]
# if in pilot mode, apply overrides according to preferences
if PILOTING:
    # force windowed mode
    if prefs.piloting['forceWindowed']:
        _fullScr = False
        # set window size
        _winSize = prefs.piloting['forcedWindowSize']

def showExpInfoDlg(expInfo):
    """
    Show participant info dialog.
    Parameters
    ==========
    expInfo : dict
        Information about this experiment.
    
    Returns
    ==========
    dict
        Information about this experiment.
    """
    # show participant info dialog
    dlg = gui.DlgFromDict(
        dictionary=expInfo, sortKeys=False, title=expName, alwaysOnTop=True
    )
    if dlg.OK == False:
        core.quit()  # user pressed cancel
    # return expInfo
    return expInfo


def setupData(expInfo, dataDir=None):
    """
    Make an ExperimentHandler to handle trials and saving.
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    dataDir : Path, str or None
        Folder to save the data to, leave as None to create a folder in the current directory.    
    Returns
    ==========
    psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    """
    # remove dialog-specific syntax from expInfo
    for key, val in expInfo.copy().items():
        newKey, _ = data.utils.parsePipeSyntax(key)
        expInfo[newKey] = expInfo.pop(key)
    
    # data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
    if dataDir is None:
        dataDir = _thisDir
    filename = u'results/%s_subj%s_order%s_%s' % (expName, expInfo['subj'], expInfo['sess_type'], expInfo['date'])
    # make sure filename is relative to dataDir
    if os.path.isabs(filename):
        dataDir = os.path.commonprefix([dataDir, filename])
        filename = os.path.relpath(filename, dataDir)
    
    # an ExperimentHandler isn't essential but helps with data saving
    thisExp = data.ExperimentHandler(
        name=expName, version='',
        extraInfo=expInfo, runtimeInfo=None,
        originPath='/Users/f0064z8/Library/CloudStorage/GoogleDrive-si2442@columbia.edu/My Drive/research/asymmetry_design/asymmetry_v9_streamline_lastrun.py',
        savePickle=True, saveWideText=True,
        dataFileName=dataDir + os.sep + filename, sortColumns='time'
    )
    thisExp.setPriority('thisRow.t', priority.CRITICAL)
    thisExp.setPriority('expName', priority.LOW)
    # return experiment handler
    return thisExp


def setupLogging(filename):
    """
    Setup a log file and tell it what level to log at.
    
    Parameters
    ==========
    filename : str or pathlib.Path
        Filename to save log file and data files as, doesn't need an extension.
    
    Returns
    ==========
    psychopy.logging.LogFile
        Text stream to receive inputs from the logging system.
    """
    # set how much information should be printed to the console / app
    if PILOTING:
        logging.console.setLevel(
            prefs.piloting['pilotConsoleLoggingLevel']
        )
    else:
        logging.console.setLevel('warning')


def setupWindow(expInfo=None, win=None):
    """
    Setup the Window
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    win : psychopy.visual.Window
        Window to setup - leave as None to create a new window.
    
    Returns
    ==========
    psychopy.visual.Window
        Window in which to run this experiment.
    """
    if PILOTING:
        logging.debug('Fullscreen settings ignored as running in pilot mode.')
    
    if win is None:
        # if not given a window to setup, make one
        win = visual.Window(
            size=_winSize, fullscr=_fullScr, screen=0,
            winType='pyglet', allowGUI=False, allowStencil=False,
            monitor='testMonitor', color=[0,0,0], colorSpace='rgb',
            backgroundImage='', backgroundFit='none',
            blendMode='avg', useFBO=True,
            units='height',
            checkTiming=False  # we're going to do this ourselves in a moment
        )
    else:
        # if we have a window, just set the attributes which are safe to set
        win.color = [0,0,0]
        win.colorSpace = 'rgb'
        win.backgroundImage = ''
        win.backgroundFit = 'none'
        win.units = 'height'
    if expInfo is not None:
        # get/measure frame rate if not already in expInfo
        if win._monitorFrameRate is None:
            win._monitorFrameRate = win.getActualFrameRate(infoMsg='Attempting to measure frame rate of screen, please wait...')
        expInfo['frameRate'] = win._monitorFrameRate
    win.hideMessage()
    # show a visual indicator if we're in piloting mode
    if PILOTING and prefs.piloting['showPilotingIndicator']:
        win.showPilotingIndicator()
    
    return win


def setupDevices(expInfo, thisExp, win):
    """
    Setup whatever devices are available (mouse, keyboard, speaker, eyetracker, etc.) and add them to 
    the device manager (deviceManager)
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    win : psychopy.visual.Window
        Window in which to run this experiment.
    Returns
    ==========
    bool
        True if completed successfully.
    """
    # --- Setup input devices ---
    ioConfig = {}
    ioSession = ioServer = eyetracker = None
    
    # store ioServer object in the device manager
    deviceManager.ioServer = ioServer
    
    # create a default keyboard (e.g. to check for escape)
    if deviceManager.getDevice('defaultKeyboard') is None:
        deviceManager.addDevice(
            deviceClass='keyboard', deviceName='defaultKeyboard', backend='ptb'
        )
    if deviceManager.getDevice('block_start_resp') is None:
        # initialise block_start_resp
        block_start_resp = deviceManager.addDevice(
            deviceClass='keyboard',
            deviceName='block_start_resp',
        )
    if deviceManager.getDevice('slider_resp') is None:
        # initialise slider_resp
        slider_resp = deviceManager.addDevice(
            deviceClass='keyboard',
            deviceName='slider_resp',
        )
    if deviceManager.getDevice('submit_resp') is None:
        # initialise submit_resp
        submit_resp = deviceManager.addDevice(
            deviceClass='keyboard',
            deviceName='submit_resp',
        )
    # return True if completed successfully
    return True

def pauseExperiment(thisExp, win=None, timers=[], playbackComponents=[]):
    """
    Pause this experiment, preventing the flow from advancing to the next routine until resumed.
    
    Parameters
    ==========
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    win : psychopy.visual.Window
        Window for this experiment.
    timers : list, tuple
        List of timers to reset once pausing is finished.
    playbackComponents : list, tuple
        List of any components with a `pause` method which need to be paused.
    """
    # if we are not paused, do nothing
    if thisExp.status != PAUSED:
        return
    
    # start a timer to figure out how long we're paused for
    pauseTimer = core.Clock()
    # pause any playback components
    for comp in playbackComponents:
        comp.pause()
    # make sure we have a keyboard
    defaultKeyboard = deviceManager.getDevice('defaultKeyboard')
    if defaultKeyboard is None:
        defaultKeyboard = deviceManager.addKeyboard(
            deviceClass='keyboard',
            deviceName='defaultKeyboard',
            backend='PsychToolbox',
        )
    # run a while loop while we wait to unpause
    while thisExp.status == PAUSED:
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=['escape']):
            endExperiment(thisExp, win=win)
        # sleep 1ms so other threads can execute
        clock.time.sleep(0.001)
    # if stop was requested while paused, quit
    if thisExp.status == FINISHED:
        endExperiment(thisExp, win=win)
    # resume any playback components
    for comp in playbackComponents:
        comp.play()
    # reset any timers
    for timer in timers:
        timer.addTime(-pauseTimer.getTime())


def run(expInfo, thisExp, win, globalClock=None, thisSession=None):
    """
    Run the experiment flow.
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    psychopy.visual.Window
        Window in which to run this experiment.
    globalClock : psychopy.core.clock.Clock or None
        Clock to get global time from - supply None to make a new one.
    thisSession : psychopy.session.Session or None
        Handle of the Session object this experiment is being run from, if any.
    """
    # mark experiment as started
    thisExp.status = STARTED
    # make sure window is set to foreground to prevent losing focus
    win.winHandle.activate()
    # make sure variables created by exec are available globally
    exec = environmenttools.setExecEnvironment(globals())
    # get device handles from dict of input devices
    ioServer = deviceManager.ioServer
    # get/create a default keyboard (e.g. to check for escape)
    defaultKeyboard = deviceManager.getDevice('defaultKeyboard')
    if defaultKeyboard is None:
        deviceManager.addDevice(
            deviceClass='keyboard', deviceName='defaultKeyboard', backend='PsychToolbox'
        )
    eyetracker = deviceManager.getDevice('eyetracker')
    # make sure we're running in the directory for this experiment
    os.chdir(_thisDir)
    # get filename from ExperimentHandler for convenience
    filename = thisExp.dataFileName
    frameTolerance = 0.001  # how close to onset before 'same' frame
    endExpNow = False  # flag for 'escape' or other condition => quit the exp
    # get frame duration from frame rate in expInfo
    if 'frameRate' in expInfo and expInfo['frameRate'] is not None:
        frameDur = 1.0 / round(expInfo['frameRate'])
    else:
        frameDur = 1.0 / 60.0  # could not measure, so guess
    
    # Start Code - component code to be run after the window creation
    
    # --- Initialize components for Routine "new_block" ---
    # Run 'Begin Experiment' code from block_start_code
    block_outcome, block_bonus = 0, 0
    #block_file = f'block_order{expInfo["block_order"]}.csv'
    sess_type_file = f'input_data/sess_type{expInfo["sess_type"]}.csv'
    block_start_text = visual.TextStim(win=win, name='block_start_text',
        text='',
        font='Arial',
        pos=(0, 0), draggable=False, height=0.05, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-1.0);
    block_start_resp = keyboard.Keyboard(deviceName='block_start_resp')
    
    # --- Initialize components for Routine "stim" ---
    targ_pos = visual.ImageStim(
        win=win,
        name='targ_pos', 
        image='default.png', mask=None, anchor='center',
        ori=0.0, pos=(0, 0), draggable=False, size=[.6,.6],
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=0.0)
    
    # --- Initialize components for Routine "delay" ---
    
    # --- Initialize components for Routine "task" ---
    # Run 'Begin Experiment' code from rand_divider_slider
    leftPressed, rightPressed, sliderMoved = .0, .0, .0
    positions = []
    debug_task_txt = ''
    img1 = visual.ImageStim(
        win=win,
        name='img1', 
        image='default.png', mask=None, anchor='center',
        ori=0.0, pos=[-.6, 0], draggable=False, size=[.5,.5],
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-1.0)
    img2 = visual.ImageStim(
        win=win,
        name='img2', 
        image='default.png', mask=None, anchor='center',
        ori=0.0, pos=[.6, 0], draggable=False, size=[.5,.5],
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-2.0)
    divider_line = visual.Line(
        win=win, name='divider_line',
        size=[.2,0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=5.0,
        colorSpace='rgb', lineColor=[1.0000, 1.0000, 1.0000], fillColor=[1.0000, 1.0000, 1.0000],
        opacity=None, depth=-3.0, interpolate=True)
    slider_line = visual.Line(
        win=win, name='slider_line',
        size=[.8,.0],
        ori=0.0, pos=[0, 0], draggable=False, anchor='center',
        lineWidth=5.0,
        colorSpace='rgb', lineColor=[1.0000, 1.0000, 1.0000], fillColor=[1.0000, 1.0000, 1.0000],
        opacity=1.0, depth=-4.0, interpolate=True)
    slider = visual.Slider(win=win, name='slider',
        startValue=None, size=[.8,.05], pos=[0,0], units=win.units,
        labels=None, ticks=(-.4,-.3,-.2,-.1,.0,.1,.2,.3,.4), granularity=0.001,
        style='slider', styleTweaks=('labels45',), opacity=1.0,
        labelColor=None, markerColor=[1.0000, -.5000, 1.0000], lineColor=None, colorSpace='rgb',
        font='Open Sans', labelHeight=0.05,
        flip=False, ori=0.0, depth=-5, readOnly=False)
    slider_resp = keyboard.Keyboard(deviceName='slider_resp')
    submit_resp = keyboard.Keyboard(deviceName='submit_resp')
    debug_task_text = visual.TextStim(win=win, name='debug_task_text',
        text='',
        font='Arial',
        pos=[0,0], draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-8.0);
    
    # --- Initialize components for Routine "feedback" ---
    # Run 'Begin Experiment' code from fb_code
    correct, outcome = .0, .0
    no_resp_text = visual.TextStim(win=win, name='no_resp_text',
        text='',
        font='Open Sans',
        pos=[0,0], draggable=False, height=0.05, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-1.0);
    Lcoin = visual.ImageStim(
        win=win,
        name='Lcoin', 
        image='default.png', mask=None, anchor='center',
        ori=0.0, pos=[-.2, 0], draggable=False, size=[.3,.3],
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-2.0)
    Rcoin = visual.ImageStim(
        win=win,
        name='Rcoin', 
        image='default.png', mask=None, anchor='center',
        ori=0.0, pos=[.2, 0], draggable=False, size=[.3,.3],
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-3.0)
    Mcoin = visual.ImageStim(
        win=win,
        name='Mcoin', 
        image='default.png', mask=None, anchor='center',
        ori=0.0, pos=(0, 0), draggable=False, size=[.3,.3],
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-4.0)
    Lcross = visual.ImageStim(
        win=win,
        name='Lcross', 
        image='default.png', mask=None, anchor='center',
        ori=0.0, pos=[-.2, 0], draggable=False, size=[.3,.3],
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-5.0)
    Rcross = visual.ImageStim(
        win=win,
        name='Rcross', 
        image='default.png', mask=None, anchor='center',
        ori=0.0, pos=[.2, 0], draggable=False, size=[.3,.3],
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-6.0)
    Mcross = visual.ImageStim(
        win=win,
        name='Mcross', 
        image='default.png', mask=None, anchor='center',
        ori=0.0, pos=(0, 0), draggable=False, size=[.3,.3],
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-7.0)
    
    # --- Initialize components for Routine "ITI" ---
    
    # --- Initialize components for Routine "end_block" ---
    block_end_text = visual.TextStim(win=win, name='block_end_text',
        text='',
        font='Arial',
        pos=(0, 0), draggable=False, height=0.05, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-1.0);
    
    # create some handy timers
    
    # global clock to track the time since experiment started
    if globalClock is None:
        # create a clock if not given one
        globalClock = core.Clock()
    if isinstance(globalClock, str):
        # if given a string, make a clock accoridng to it
        if globalClock == 'float':
            # get timestamps as a simple value
            globalClock = core.Clock(format='float')
        elif globalClock == 'iso':
            # get timestamps in ISO format
            globalClock = core.Clock(format='%Y-%m-%d_%H:%M:%S.%f%z')
        else:
            # get timestamps in a custom format
            globalClock = core.Clock(format=globalClock)
    if ioServer is not None:
        ioServer.syncClock(globalClock)
    logging.setDefaultClock(globalClock)
    # routine timer to track time remaining of each (possibly non-slip) routine
    routineTimer = core.Clock()
    win.flip()  # flip window to reset last flip timer
    # store the exact time the global clock started
    expInfo['expStart'] = data.getDateStr(
        format='%Y-%m-%d %Hh%M.%S.%f %z', fractionalSecondDigits=6
    )
    
    # set up handler to look after randomisation of conditions etc
    blocks = data.TrialHandler2(
        name='blocks',
        nReps=1.0, 
        method='sequential', 
        extraInfo=expInfo, 
        originPath=-1, 
        trialList=data.importConditions(sess_type_file), 
        seed=None, 
    )
    thisExp.addLoop(blocks)  # add the loop to the experiment
    thisBlock = blocks.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisBlock.rgb)
    if thisBlock != None:
        for paramName in thisBlock:
            globals()[paramName] = thisBlock[paramName]
    
    for thisBlock in blocks:
        currentLoop = blocks
        thisExp.timestampOnFlip(win, 'thisRow.t', format=globalClock.format)
        # abbreviate parameter names if possible (e.g. rgb = thisBlock.rgb)
        if thisBlock != None:
            for paramName in thisBlock:
                globals()[paramName] = thisBlock[paramName]
        
        # --- Prepare to start Routine "new_block" ---
        # create an object to store info about Routine new_block
        new_block = data.Routine(
            name='new_block',
            components=[block_start_text, block_start_resp],
        )
        new_block.status = NOT_STARTED
        continueRoutine = True
        # update component parameters for each repeat
        # Run 'Begin Routine' code from block_start_code
        block_outcome, block_bonus = 0, 0
        block_start_txt = f'Press space to begin block {blockN}'
        
        trial_rows = [i + (blockN-1) * 40 for i in range(40)]
        block_start_text.setText(block_start_txt)
        # create starting attributes for block_start_resp
        block_start_resp.keys = []
        block_start_resp.rt = []
        _block_start_resp_allKeys = []
        # store start times for new_block
        new_block.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
        new_block.tStart = globalClock.getTime(format='float')
        new_block.status = STARTED
        thisExp.addData('new_block.started', new_block.tStart)
        new_block.maxDuration = None
        # keep track of which components have finished
        new_blockComponents = new_block.components
        for thisComponent in new_block.components:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "new_block" ---
        # if trial has changed, end Routine now
        if isinstance(blocks, data.TrialHandler2) and thisBlock.thisN != blocks.thisTrial.thisN:
            continueRoutine = False
        new_block.forceEnded = routineForceEnded = not continueRoutine
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *block_start_text* updates
            
            # if block_start_text is starting this frame...
            if block_start_text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                block_start_text.frameNStart = frameN  # exact frame index
                block_start_text.tStart = t  # local t and not account for scr refresh
                block_start_text.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(block_start_text, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'block_start_text.started')
                # update status
                block_start_text.status = STARTED
                block_start_text.setAutoDraw(True)
            
            # if block_start_text is active this frame...
            if block_start_text.status == STARTED:
                # update params
                pass
            
            # *block_start_resp* updates
            waitOnFlip = False
            
            # if block_start_resp is starting this frame...
            if block_start_resp.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                block_start_resp.frameNStart = frameN  # exact frame index
                block_start_resp.tStart = t  # local t and not account for scr refresh
                block_start_resp.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(block_start_resp, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'block_start_resp.started')
                # update status
                block_start_resp.status = STARTED
                # keyboard checking is just starting
                waitOnFlip = True
                win.callOnFlip(block_start_resp.clock.reset)  # t=0 on next screen flip
                win.callOnFlip(block_start_resp.clearEvents, eventType='keyboard')  # clear events on next screen flip
            if block_start_resp.status == STARTED and not waitOnFlip:
                theseKeys = block_start_resp.getKeys(keyList=['space'], ignoreKeys=["escape"], waitRelease=False)
                _block_start_resp_allKeys.extend(theseKeys)
                if len(_block_start_resp_allKeys):
                    block_start_resp.keys = _block_start_resp_allKeys[-1].name  # just the last key pressed
                    block_start_resp.rt = _block_start_resp_allKeys[-1].rt
                    block_start_resp.duration = _block_start_resp_allKeys[-1].duration
                    # a response ends the routine
                    continueRoutine = False
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, win=win)
                return
            # pause experiment here if requested
            if thisExp.status == PAUSED:
                pauseExperiment(
                    thisExp=thisExp, 
                    win=win, 
                    timers=[routineTimer], 
                    playbackComponents=[]
                )
                # skip the frame we paused on
                continue
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                new_block.forceEnded = routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in new_block.components:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "new_block" ---
        for thisComponent in new_block.components:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # store stop times for new_block
        new_block.tStop = globalClock.getTime(format='float')
        new_block.tStopRefresh = tThisFlipGlobal
        thisExp.addData('new_block.stopped', new_block.tStop)
        # check responses
        if block_start_resp.keys in ['', [], None]:  # No response was made
            block_start_resp.keys = None
        blocks.addData('block_start_resp.keys',block_start_resp.keys)
        if block_start_resp.keys != None:  # we had a response
            blocks.addData('block_start_resp.rt', block_start_resp.rt)
            blocks.addData('block_start_resp.duration', block_start_resp.duration)
        # the Routine "new_block" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # set up handler to look after randomisation of conditions etc
        trials = data.TrialHandler2(
            name='trials',
            nReps=1.0, 
            method='random', 
            extraInfo=expInfo, 
            originPath=-1, 
            trialList=data.importConditions(
            'input_data/trials.csv', 
            selection=trial_rows
        )
        , 
            seed=None, 
        )
        thisExp.addLoop(trials)  # add the loop to the experiment
        thisTrial = trials.trialList[0]  # so we can initialise stimuli with some values
        # abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
        if thisTrial != None:
            for paramName in thisTrial:
                globals()[paramName] = thisTrial[paramName]
        if thisSession is not None:
            # if running in a Session with a Liaison client, send data up to now
            thisSession.sendExperimentData()
        
        for thisTrial in trials:
            currentLoop = trials
            thisExp.timestampOnFlip(win, 'thisRow.t', format=globalClock.format)
            if thisSession is not None:
                # if running in a Session with a Liaison client, send data up to now
                thisSession.sendExperimentData()
            # abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
            if thisTrial != None:
                for paramName in thisTrial:
                    globals()[paramName] = thisTrial[paramName]
            
            # --- Prepare to start Routine "stim" ---
            # create an object to store info about Routine stim
            stim = data.Routine(
                name='stim',
                components=[targ_pos],
            )
            stim.status = NOT_STARTED
            continueRoutine = True
            # update component parameters for each repeat
            targ_pos.setImage(target_file)
            # store start times for stim
            stim.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
            stim.tStart = globalClock.getTime(format='float')
            stim.status = STARTED
            thisExp.addData('stim.started', stim.tStart)
            stim.maxDuration = None
            # keep track of which components have finished
            stimComponents = stim.components
            for thisComponent in stim.components:
                thisComponent.tStart = None
                thisComponent.tStop = None
                thisComponent.tStartRefresh = None
                thisComponent.tStopRefresh = None
                if hasattr(thisComponent, 'status'):
                    thisComponent.status = NOT_STARTED
            # reset timers
            t = 0
            _timeToFirstFrame = win.getFutureFlipTime(clock="now")
            frameN = -1
            
            # --- Run Routine "stim" ---
            # if trial has changed, end Routine now
            if isinstance(trials, data.TrialHandler2) and thisTrial.thisN != trials.thisTrial.thisN:
                continueRoutine = False
            stim.forceEnded = routineForceEnded = not continueRoutine
            while continueRoutine and routineTimer.getTime() < 1.0:
                # get current time
                t = routineTimer.getTime()
                tThisFlip = win.getFutureFlipTime(clock=routineTimer)
                tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                # update/draw components on each frame
                
                # *targ_pos* updates
                
                # if targ_pos is starting this frame...
                if targ_pos.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    targ_pos.frameNStart = frameN  # exact frame index
                    targ_pos.tStart = t  # local t and not account for scr refresh
                    targ_pos.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(targ_pos, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'targ_pos.started')
                    # update status
                    targ_pos.status = STARTED
                    targ_pos.setAutoDraw(True)
                
                # if targ_pos is active this frame...
                if targ_pos.status == STARTED:
                    # update params
                    pass
                
                # if targ_pos is stopping this frame...
                if targ_pos.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > targ_pos.tStartRefresh + 1-frameTolerance:
                        # keep track of stop time/frame for later
                        targ_pos.tStop = t  # not accounting for scr refresh
                        targ_pos.tStopRefresh = tThisFlipGlobal  # on global time
                        targ_pos.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'targ_pos.stopped')
                        # update status
                        targ_pos.status = FINISHED
                        targ_pos.setAutoDraw(False)
                
                # check for quit (typically the Esc key)
                if defaultKeyboard.getKeys(keyList=["escape"]):
                    thisExp.status = FINISHED
                if thisExp.status == FINISHED or endExpNow:
                    endExperiment(thisExp, win=win)
                    return
                # pause experiment here if requested
                if thisExp.status == PAUSED:
                    pauseExperiment(
                        thisExp=thisExp, 
                        win=win, 
                        timers=[routineTimer], 
                        playbackComponents=[]
                    )
                    # skip the frame we paused on
                    continue
                
                # check if all components have finished
                if not continueRoutine:  # a component has requested a forced-end of Routine
                    stim.forceEnded = routineForceEnded = True
                    break
                continueRoutine = False  # will revert to True if at least one component still running
                for thisComponent in stim.components:
                    if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                        continueRoutine = True
                        break  # at least one component has not yet finished
                
                # refresh the screen
                if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                    win.flip()
            
            # --- Ending Routine "stim" ---
            for thisComponent in stim.components:
                if hasattr(thisComponent, "setAutoDraw"):
                    thisComponent.setAutoDraw(False)
            # store stop times for stim
            stim.tStop = globalClock.getTime(format='float')
            stim.tStopRefresh = tThisFlipGlobal
            thisExp.addData('stim.stopped', stim.tStop)
            # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
            if stim.maxDurationReached:
                routineTimer.addTime(-stim.maxDuration)
            elif stim.forceEnded:
                routineTimer.reset()
            else:
                routineTimer.addTime(-1.000000)
            
            # --- Prepare to start Routine "delay" ---
            # create an object to store info about Routine delay
            delay = data.Routine(
                name='delay',
                components=[],
            )
            delay.status = NOT_STARTED
            continueRoutine = True
            # update component parameters for each repeat
            # store start times for delay
            delay.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
            delay.tStart = globalClock.getTime(format='float')
            delay.status = STARTED
            thisExp.addData('delay.started', delay.tStart)
            delay.maxDuration = 1.5
            # keep track of which components have finished
            delayComponents = delay.components
            for thisComponent in delay.components:
                thisComponent.tStart = None
                thisComponent.tStop = None
                thisComponent.tStartRefresh = None
                thisComponent.tStopRefresh = None
                if hasattr(thisComponent, 'status'):
                    thisComponent.status = NOT_STARTED
            # reset timers
            t = 0
            _timeToFirstFrame = win.getFutureFlipTime(clock="now")
            frameN = -1
            
            # --- Run Routine "delay" ---
            # if trial has changed, end Routine now
            if isinstance(trials, data.TrialHandler2) and thisTrial.thisN != trials.thisTrial.thisN:
                continueRoutine = False
            delay.forceEnded = routineForceEnded = not continueRoutine
            while continueRoutine and routineTimer.getTime() < 1.5:
                # get current time
                t = routineTimer.getTime()
                tThisFlip = win.getFutureFlipTime(clock=routineTimer)
                tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                # update/draw components on each frame
                # is it time to end the Routine? (based on local clock)
                if tThisFlip > delay.maxDuration-frameTolerance:
                    delay.maxDurationReached = True
                    continueRoutine = False
                
                # check for quit (typically the Esc key)
                if defaultKeyboard.getKeys(keyList=["escape"]):
                    thisExp.status = FINISHED
                if thisExp.status == FINISHED or endExpNow:
                    endExperiment(thisExp, win=win)
                    return
                # pause experiment here if requested
                if thisExp.status == PAUSED:
                    pauseExperiment(
                        thisExp=thisExp, 
                        win=win, 
                        timers=[routineTimer], 
                        playbackComponents=[]
                    )
                    # skip the frame we paused on
                    continue
                
                # check if all components have finished
                if not continueRoutine:  # a component has requested a forced-end of Routine
                    delay.forceEnded = routineForceEnded = True
                    break
                continueRoutine = False  # will revert to True if at least one component still running
                for thisComponent in delay.components:
                    if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                        continueRoutine = True
                        break  # at least one component has not yet finished
                
                # refresh the screen
                if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                    win.flip()
            
            # --- Ending Routine "delay" ---
            for thisComponent in delay.components:
                if hasattr(thisComponent, "setAutoDraw"):
                    thisComponent.setAutoDraw(False)
            # store stop times for delay
            delay.tStop = globalClock.getTime(format='float')
            delay.tStopRefresh = tThisFlipGlobal
            thisExp.addData('delay.stopped', delay.tStop)
            # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
            if delay.maxDurationReached:
                routineTimer.addTime(-delay.maxDuration)
            elif delay.forceEnded:
                routineTimer.reset()
            else:
                routineTimer.addTime(-1.500000)
            
            # --- Prepare to start Routine "task" ---
            # create an object to store info about Routine task
            task = data.Routine(
                name='task',
                components=[img1, img2, divider_line, slider_line, slider, slider_resp, submit_resp, debug_task_text],
            )
            task.status = NOT_STARTED
            continueRoutine = True
            # update component parameters for each repeat
            # Run 'Begin Routine' code from rand_divider_slider
            leftPressed, rightPressed, sliderMoved = .0, .0, .0
            positions = []
            debug_task_txt = f'Trial {trial_key+1}'
            
            # to collect responses
            kb.clearEvents()
            img1.setImage(img1_file)
            img2.setImage(img2_file)
            divider_line.setPos([disp_div, 0])
            slider.reset()
            # create starting attributes for slider_resp
            slider_resp.keys = []
            slider_resp.rt = []
            _slider_resp_allKeys = []
            # create starting attributes for submit_resp
            submit_resp.keys = []
            submit_resp.rt = []
            _submit_resp_allKeys = []
            debug_task_text.setText(debug_task_txt)
            # store start times for task
            task.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
            task.tStart = globalClock.getTime(format='float')
            task.status = STARTED
            thisExp.addData('task.started', task.tStart)
            task.maxDuration = None
            # keep track of which components have finished
            taskComponents = task.components
            for thisComponent in task.components:
                thisComponent.tStart = None
                thisComponent.tStop = None
                thisComponent.tStartRefresh = None
                thisComponent.tStopRefresh = None
                if hasattr(thisComponent, 'status'):
                    thisComponent.status = NOT_STARTED
            # reset timers
            t = 0
            _timeToFirstFrame = win.getFutureFlipTime(clock="now")
            frameN = -1
            
            # --- Run Routine "task" ---
            # if trial has changed, end Routine now
            if isinstance(trials, data.TrialHandler2) and thisTrial.thisN != trials.thisTrial.thisN:
                continueRoutine = False
            task.forceEnded = routineForceEnded = not continueRoutine
            while continueRoutine and routineTimer.getTime() < 3.0:
                # get current time
                t = routineTimer.getTime()
                tThisFlip = win.getFutureFlipTime(clock=routineTimer)
                tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                # update/draw components on each frame
                # Run 'Each Frame' code from rand_divider_slider
                # for some reason psychopy doesnt remember this assignment
                if slider.markerPos == None:
                    slider.markerPos = anchor
                
                if (not leftPressed and kb.getKeys(['left'], waitRelease=False, clear=False))\
                or (leftPressed and not kb.getKeys(['left'], waitRelease=True, clear=False)):
                    leftPressed, sliderMoved = 1, 1
                    if slider.markerPos >= (-.4 + .01):
                        slider.markerPos -= 0.01
                elif not leftPressed or kb.getKeys(['left'], waitRelease=True, clear=True):
                    leftPressed = 0
                    
                if (not rightPressed and kb.getKeys(['right'], waitRelease=False, clear=False))\
                or (rightPressed and not kb.getKeys(['right'], waitRelease=True, clear=False)):
                    rightPressed, sliderMoved = 1, 1    
                    if slider.markerPos <= (.4 - .01):
                        slider.markerPos += 0.01
                elif not rightPressed or kb.getKeys(['right'], waitRelease=True, clear=True):
                    rightPressed = 0
                
                positions.append(slider.markerPos)
                
                # check psychopy_notes for link
                
                ################################## debugging
                #if (slider.markerPos >= disp_divider and target >= disp_divider)\
                #or (slider.markerPos <= disp_divider and target <= disp_divider):
                #    correct = 1
                #    outcome = 2 if context == 'rew' else 1
                #
                #    error = abs(target - slider.markerPos)
                #    addn_pts = 1 - 1.5 * error
                #    
                #else:
                #    correct = -1
                #    outcome = -1 if context == 'rew' else -2
                #    error = abs(target - slider.markerPos)
                #    addn_pts = - 1.5 * error 
                #    
                #net_reward = outcome + addn_pts
                #
                #debug_task_txt  = f'{addn_pts:.03f}'
                
                ################################## end of debugging
                
                # *img1* updates
                
                # if img1 is starting this frame...
                if img1.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    img1.frameNStart = frameN  # exact frame index
                    img1.tStart = t  # local t and not account for scr refresh
                    img1.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(img1, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'img1.started')
                    # update status
                    img1.status = STARTED
                    img1.setAutoDraw(True)
                
                # if img1 is active this frame...
                if img1.status == STARTED:
                    # update params
                    pass
                
                # if img1 is stopping this frame...
                if img1.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > img1.tStartRefresh + 3-frameTolerance:
                        # keep track of stop time/frame for later
                        img1.tStop = t  # not accounting for scr refresh
                        img1.tStopRefresh = tThisFlipGlobal  # on global time
                        img1.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'img1.stopped')
                        # update status
                        img1.status = FINISHED
                        img1.setAutoDraw(False)
                
                # *img2* updates
                
                # if img2 is starting this frame...
                if img2.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    img2.frameNStart = frameN  # exact frame index
                    img2.tStart = t  # local t and not account for scr refresh
                    img2.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(img2, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'img2.started')
                    # update status
                    img2.status = STARTED
                    img2.setAutoDraw(True)
                
                # if img2 is active this frame...
                if img2.status == STARTED:
                    # update params
                    pass
                
                # if img2 is stopping this frame...
                if img2.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > img2.tStartRefresh + 3-frameTolerance:
                        # keep track of stop time/frame for later
                        img2.tStop = t  # not accounting for scr refresh
                        img2.tStopRefresh = tThisFlipGlobal  # on global time
                        img2.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'img2.stopped')
                        # update status
                        img2.status = FINISHED
                        img2.setAutoDraw(False)
                
                # *divider_line* updates
                
                # if divider_line is starting this frame...
                if divider_line.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    divider_line.frameNStart = frameN  # exact frame index
                    divider_line.tStart = t  # local t and not account for scr refresh
                    divider_line.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(divider_line, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'divider_line.started')
                    # update status
                    divider_line.status = STARTED
                    divider_line.setAutoDraw(True)
                
                # if divider_line is active this frame...
                if divider_line.status == STARTED:
                    # update params
                    pass
                
                # if divider_line is stopping this frame...
                if divider_line.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > divider_line.tStartRefresh + 3-frameTolerance:
                        # keep track of stop time/frame for later
                        divider_line.tStop = t  # not accounting for scr refresh
                        divider_line.tStopRefresh = tThisFlipGlobal  # on global time
                        divider_line.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'divider_line.stopped')
                        # update status
                        divider_line.status = FINISHED
                        divider_line.setAutoDraw(False)
                
                # *slider_line* updates
                
                # if slider_line is starting this frame...
                if slider_line.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    slider_line.frameNStart = frameN  # exact frame index
                    slider_line.tStart = t  # local t and not account for scr refresh
                    slider_line.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(slider_line, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'slider_line.started')
                    # update status
                    slider_line.status = STARTED
                    slider_line.setAutoDraw(True)
                
                # if slider_line is active this frame...
                if slider_line.status == STARTED:
                    # update params
                    pass
                
                # if slider_line is stopping this frame...
                if slider_line.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > slider_line.tStartRefresh + 3-frameTolerance:
                        # keep track of stop time/frame for later
                        slider_line.tStop = t  # not accounting for scr refresh
                        slider_line.tStopRefresh = tThisFlipGlobal  # on global time
                        slider_line.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'slider_line.stopped')
                        # update status
                        slider_line.status = FINISHED
                        slider_line.setAutoDraw(False)
                
                # *slider* updates
                
                # if slider is starting this frame...
                if slider.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    slider.frameNStart = frameN  # exact frame index
                    slider.tStart = t  # local t and not account for scr refresh
                    slider.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(slider, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'slider.started')
                    # update status
                    slider.status = STARTED
                    slider.setAutoDraw(True)
                
                # if slider is active this frame...
                if slider.status == STARTED:
                    # update params
                    pass
                
                # if slider is stopping this frame...
                if slider.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > slider.tStartRefresh + 3-frameTolerance:
                        # keep track of stop time/frame for later
                        slider.tStop = t  # not accounting for scr refresh
                        slider.tStopRefresh = tThisFlipGlobal  # on global time
                        slider.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'slider.stopped')
                        # update status
                        slider.status = FINISHED
                        slider.setAutoDraw(False)
                
                # *slider_resp* updates
                waitOnFlip = False
                
                # if slider_resp is starting this frame...
                if slider_resp.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    slider_resp.frameNStart = frameN  # exact frame index
                    slider_resp.tStart = t  # local t and not account for scr refresh
                    slider_resp.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(slider_resp, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'slider_resp.started')
                    # update status
                    slider_resp.status = STARTED
                    # keyboard checking is just starting
                    waitOnFlip = True
                    win.callOnFlip(slider_resp.clock.reset)  # t=0 on next screen flip
                
                # if slider_resp is stopping this frame...
                if slider_resp.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > slider_resp.tStartRefresh + 3-frameTolerance:
                        # keep track of stop time/frame for later
                        slider_resp.tStop = t  # not accounting for scr refresh
                        slider_resp.tStopRefresh = tThisFlipGlobal  # on global time
                        slider_resp.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'slider_resp.stopped')
                        # update status
                        slider_resp.status = FINISHED
                        slider_resp.status = FINISHED
                if slider_resp.status == STARTED and not waitOnFlip:
                    theseKeys = slider_resp.getKeys(keyList=['left','right'], ignoreKeys=["escape"], waitRelease=False)
                    _slider_resp_allKeys.extend(theseKeys)
                    if len(_slider_resp_allKeys):
                        slider_resp.keys = [key.name for key in _slider_resp_allKeys]  # storing all keys
                        slider_resp.rt = [key.rt for key in _slider_resp_allKeys]
                        slider_resp.duration = [key.duration for key in _slider_resp_allKeys]
                
                # *submit_resp* updates
                waitOnFlip = False
                
                # if submit_resp is starting this frame...
                if submit_resp.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    submit_resp.frameNStart = frameN  # exact frame index
                    submit_resp.tStart = t  # local t and not account for scr refresh
                    submit_resp.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(submit_resp, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'submit_resp.started')
                    # update status
                    submit_resp.status = STARTED
                    # keyboard checking is just starting
                    waitOnFlip = True
                    win.callOnFlip(submit_resp.clock.reset)  # t=0 on next screen flip
                    win.callOnFlip(submit_resp.clearEvents, eventType='keyboard')  # clear events on next screen flip
                
                # if submit_resp is stopping this frame...
                if submit_resp.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > submit_resp.tStartRefresh + 3-frameTolerance:
                        # keep track of stop time/frame for later
                        submit_resp.tStop = t  # not accounting for scr refresh
                        submit_resp.tStopRefresh = tThisFlipGlobal  # on global time
                        submit_resp.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'submit_resp.stopped')
                        # update status
                        submit_resp.status = FINISHED
                        submit_resp.status = FINISHED
                if submit_resp.status == STARTED and not waitOnFlip:
                    theseKeys = submit_resp.getKeys(keyList=['space'], ignoreKeys=["escape"], waitRelease=False)
                    _submit_resp_allKeys.extend(theseKeys)
                    if len(_submit_resp_allKeys):
                        submit_resp.keys = _submit_resp_allKeys[-1].name  # just the last key pressed
                        submit_resp.rt = _submit_resp_allKeys[-1].rt
                        submit_resp.duration = _submit_resp_allKeys[-1].duration
                        # a response ends the routine
                        continueRoutine = False
                
                # *debug_task_text* updates
                
                # if debug_task_text is starting this frame...
                if debug_task_text.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    debug_task_text.frameNStart = frameN  # exact frame index
                    debug_task_text.tStart = t  # local t and not account for scr refresh
                    debug_task_text.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(debug_task_text, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'debug_task_text.started')
                    # update status
                    debug_task_text.status = STARTED
                    debug_task_text.setAutoDraw(True)
                
                # if debug_task_text is active this frame...
                if debug_task_text.status == STARTED:
                    # update params
                    debug_task_text.setPos([0,-.2], log=False)
                
                # if debug_task_text is stopping this frame...
                if debug_task_text.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > debug_task_text.tStartRefresh + 3-frameTolerance:
                        # keep track of stop time/frame for later
                        debug_task_text.tStop = t  # not accounting for scr refresh
                        debug_task_text.tStopRefresh = tThisFlipGlobal  # on global time
                        debug_task_text.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'debug_task_text.stopped')
                        # update status
                        debug_task_text.status = FINISHED
                        debug_task_text.setAutoDraw(False)
                
                # check for quit (typically the Esc key)
                if defaultKeyboard.getKeys(keyList=["escape"]):
                    thisExp.status = FINISHED
                if thisExp.status == FINISHED or endExpNow:
                    endExperiment(thisExp, win=win)
                    return
                # pause experiment here if requested
                if thisExp.status == PAUSED:
                    pauseExperiment(
                        thisExp=thisExp, 
                        win=win, 
                        timers=[routineTimer], 
                        playbackComponents=[]
                    )
                    # skip the frame we paused on
                    continue
                
                # check if all components have finished
                if not continueRoutine:  # a component has requested a forced-end of Routine
                    task.forceEnded = routineForceEnded = True
                    break
                continueRoutine = False  # will revert to True if at least one component still running
                for thisComponent in task.components:
                    if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                        continueRoutine = True
                        break  # at least one component has not yet finished
                
                # refresh the screen
                if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                    win.flip()
            
            # --- Ending Routine "task" ---
            for thisComponent in task.components:
                if hasattr(thisComponent, "setAutoDraw"):
                    thisComponent.setAutoDraw(False)
            # store stop times for task
            task.tStop = globalClock.getTime(format='float')
            task.tStopRefresh = tThisFlipGlobal
            thisExp.addData('task.stopped', task.tStop)
            # Run 'End Routine' code from rand_divider_slider
            thisExp.addData('positions', positions)
            trials.addData('slider.response', slider.getRating())
            trials.addData('slider.rt', slider.getRT())
            # check responses
            if slider_resp.keys in ['', [], None]:  # No response was made
                slider_resp.keys = None
            trials.addData('slider_resp.keys',slider_resp.keys)
            if slider_resp.keys != None:  # we had a response
                trials.addData('slider_resp.rt', slider_resp.rt)
                trials.addData('slider_resp.duration', slider_resp.duration)
            # check responses
            if submit_resp.keys in ['', [], None]:  # No response was made
                submit_resp.keys = None
            trials.addData('submit_resp.keys',submit_resp.keys)
            if submit_resp.keys != None:  # we had a response
                trials.addData('submit_resp.rt', submit_resp.rt)
                trials.addData('submit_resp.duration', submit_resp.duration)
            # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
            if task.maxDurationReached:
                routineTimer.addTime(-task.maxDuration)
            elif task.forceEnded:
                routineTimer.reset()
            else:
                routineTimer.addTime(-3.000000)
            
            # --- Prepare to start Routine "feedback" ---
            # create an object to store info about Routine feedback
            feedback = data.Routine(
                name='feedback',
                components=[no_resp_text, Lcoin, Rcoin, Mcoin, Lcross, Rcross, Mcross],
            )
            feedback.status = NOT_STARTED
            continueRoutine = True
            # update component parameters for each repeat
            # Run 'Begin Routine' code from fb_code
            # C0F1 = Curve penalty, Flat reward
            # stimVal comes from file which depends on experimenter input
            valence = val_C0F1 if stimVal == 'val_C0F1' else val_C1F0
                
            no_resp_txt = ''
            coinLR, coinM, crossLR, crossM = 0, 0, 0, 0
            correct, outcome = 0, 0
            
            if sliderMoved and not submit_resp.keys == None and 'space' in submit_resp.keys:
                
                if (slider.markerPos >= disp_div and target_pos >= disp_div)\
                or (slider.markerPos <= disp_div and target_pos <= disp_div):
                    
                    correct = 1
                    outcome = 2 if valence == 'rew' else 1
                    if outcome == 2:
                        coinLR = 1
                    elif outcome == 1:
                        coinM = 1
                        
                else:
                    
                    correct = -1
                    outcome = -1 if valence == 'rew' else -2
                    if outcome == -2:
                        coinLR, crossLR = 1, 1
                    elif outcome == -1:
                        coinM, crossM = 1, 1   
                
                block_outcome += outcome
                bonus = 1 - (target_pos - slider.markerPos)**2
                block_bonus += bonus
                
            elif not sliderMoved:
                no_resp_txt = 'You must move the slider!'
                
            else:
                no_resp_txt = 'Respond faster!'
            
            no_resp_text.setText(no_resp_txt)
            Lcoin.setImage('input_data/stims/coin.png')
            Rcoin.setImage('input_data/stims/coin.png')
            Mcoin.setImage('input_data/stims/coin.png')
            Lcross.setImage('input_data/stims/cross.png')
            Rcross.setImage('input_data/stims/cross.png')
            Mcross.setImage('input_data/stims/cross.png')
            # store start times for feedback
            feedback.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
            feedback.tStart = globalClock.getTime(format='float')
            feedback.status = STARTED
            thisExp.addData('feedback.started', feedback.tStart)
            feedback.maxDuration = 1
            # keep track of which components have finished
            feedbackComponents = feedback.components
            for thisComponent in feedback.components:
                thisComponent.tStart = None
                thisComponent.tStop = None
                thisComponent.tStartRefresh = None
                thisComponent.tStopRefresh = None
                if hasattr(thisComponent, 'status'):
                    thisComponent.status = NOT_STARTED
            # reset timers
            t = 0
            _timeToFirstFrame = win.getFutureFlipTime(clock="now")
            frameN = -1
            
            # --- Run Routine "feedback" ---
            # if trial has changed, end Routine now
            if isinstance(trials, data.TrialHandler2) and thisTrial.thisN != trials.thisTrial.thisN:
                continueRoutine = False
            feedback.forceEnded = routineForceEnded = not continueRoutine
            while continueRoutine:
                # get current time
                t = routineTimer.getTime()
                tThisFlip = win.getFutureFlipTime(clock=routineTimer)
                tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                # update/draw components on each frame
                # is it time to end the Routine? (based on local clock)
                if tThisFlip > feedback.maxDuration-frameTolerance:
                    feedback.maxDurationReached = True
                    continueRoutine = False
                
                # *no_resp_text* updates
                
                # if no_resp_text is starting this frame...
                if no_resp_text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    no_resp_text.frameNStart = frameN  # exact frame index
                    no_resp_text.tStart = t  # local t and not account for scr refresh
                    no_resp_text.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(no_resp_text, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'no_resp_text.started')
                    # update status
                    no_resp_text.status = STARTED
                    no_resp_text.setAutoDraw(True)
                
                # if no_resp_text is active this frame...
                if no_resp_text.status == STARTED:
                    # update params
                    pass
                
                # if no_resp_text is stopping this frame...
                if no_resp_text.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > no_resp_text.tStartRefresh + 1-frameTolerance:
                        # keep track of stop time/frame for later
                        no_resp_text.tStop = t  # not accounting for scr refresh
                        no_resp_text.tStopRefresh = tThisFlipGlobal  # on global time
                        no_resp_text.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'no_resp_text.stopped')
                        # update status
                        no_resp_text.status = FINISHED
                        no_resp_text.setAutoDraw(False)
                
                # *Lcoin* updates
                
                # if Lcoin is starting this frame...
                if Lcoin.status == NOT_STARTED and coinLR:
                    # keep track of start time/frame for later
                    Lcoin.frameNStart = frameN  # exact frame index
                    Lcoin.tStart = t  # local t and not account for scr refresh
                    Lcoin.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(Lcoin, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'Lcoin.started')
                    # update status
                    Lcoin.status = STARTED
                    Lcoin.setAutoDraw(True)
                
                # if Lcoin is active this frame...
                if Lcoin.status == STARTED:
                    # update params
                    pass
                
                # if Lcoin is stopping this frame...
                if Lcoin.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > Lcoin.tStartRefresh + 1-frameTolerance:
                        # keep track of stop time/frame for later
                        Lcoin.tStop = t  # not accounting for scr refresh
                        Lcoin.tStopRefresh = tThisFlipGlobal  # on global time
                        Lcoin.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'Lcoin.stopped')
                        # update status
                        Lcoin.status = FINISHED
                        Lcoin.setAutoDraw(False)
                
                # *Rcoin* updates
                
                # if Rcoin is starting this frame...
                if Rcoin.status == NOT_STARTED and coinLR:
                    # keep track of start time/frame for later
                    Rcoin.frameNStart = frameN  # exact frame index
                    Rcoin.tStart = t  # local t and not account for scr refresh
                    Rcoin.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(Rcoin, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'Rcoin.started')
                    # update status
                    Rcoin.status = STARTED
                    Rcoin.setAutoDraw(True)
                
                # if Rcoin is active this frame...
                if Rcoin.status == STARTED:
                    # update params
                    pass
                
                # if Rcoin is stopping this frame...
                if Rcoin.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > Rcoin.tStartRefresh + 1-frameTolerance:
                        # keep track of stop time/frame for later
                        Rcoin.tStop = t  # not accounting for scr refresh
                        Rcoin.tStopRefresh = tThisFlipGlobal  # on global time
                        Rcoin.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'Rcoin.stopped')
                        # update status
                        Rcoin.status = FINISHED
                        Rcoin.setAutoDraw(False)
                
                # *Mcoin* updates
                
                # if Mcoin is starting this frame...
                if Mcoin.status == NOT_STARTED and coinM:
                    # keep track of start time/frame for later
                    Mcoin.frameNStart = frameN  # exact frame index
                    Mcoin.tStart = t  # local t and not account for scr refresh
                    Mcoin.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(Mcoin, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'Mcoin.started')
                    # update status
                    Mcoin.status = STARTED
                    Mcoin.setAutoDraw(True)
                
                # if Mcoin is active this frame...
                if Mcoin.status == STARTED:
                    # update params
                    pass
                
                # if Mcoin is stopping this frame...
                if Mcoin.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > Mcoin.tStartRefresh + 1-frameTolerance:
                        # keep track of stop time/frame for later
                        Mcoin.tStop = t  # not accounting for scr refresh
                        Mcoin.tStopRefresh = tThisFlipGlobal  # on global time
                        Mcoin.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'Mcoin.stopped')
                        # update status
                        Mcoin.status = FINISHED
                        Mcoin.setAutoDraw(False)
                
                # *Lcross* updates
                
                # if Lcross is starting this frame...
                if Lcross.status == NOT_STARTED and crossLR:
                    # keep track of start time/frame for later
                    Lcross.frameNStart = frameN  # exact frame index
                    Lcross.tStart = t  # local t and not account for scr refresh
                    Lcross.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(Lcross, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'Lcross.started')
                    # update status
                    Lcross.status = STARTED
                    Lcross.setAutoDraw(True)
                
                # if Lcross is active this frame...
                if Lcross.status == STARTED:
                    # update params
                    pass
                
                # if Lcross is stopping this frame...
                if Lcross.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > Lcross.tStartRefresh + 1-frameTolerance:
                        # keep track of stop time/frame for later
                        Lcross.tStop = t  # not accounting for scr refresh
                        Lcross.tStopRefresh = tThisFlipGlobal  # on global time
                        Lcross.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'Lcross.stopped')
                        # update status
                        Lcross.status = FINISHED
                        Lcross.setAutoDraw(False)
                
                # *Rcross* updates
                
                # if Rcross is starting this frame...
                if Rcross.status == NOT_STARTED and crossLR:
                    # keep track of start time/frame for later
                    Rcross.frameNStart = frameN  # exact frame index
                    Rcross.tStart = t  # local t and not account for scr refresh
                    Rcross.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(Rcross, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'Rcross.started')
                    # update status
                    Rcross.status = STARTED
                    Rcross.setAutoDraw(True)
                
                # if Rcross is active this frame...
                if Rcross.status == STARTED:
                    # update params
                    pass
                
                # if Rcross is stopping this frame...
                if Rcross.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > Rcross.tStartRefresh + 1-frameTolerance:
                        # keep track of stop time/frame for later
                        Rcross.tStop = t  # not accounting for scr refresh
                        Rcross.tStopRefresh = tThisFlipGlobal  # on global time
                        Rcross.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'Rcross.stopped')
                        # update status
                        Rcross.status = FINISHED
                        Rcross.setAutoDraw(False)
                
                # *Mcross* updates
                
                # if Mcross is starting this frame...
                if Mcross.status == NOT_STARTED and crossM:
                    # keep track of start time/frame for later
                    Mcross.frameNStart = frameN  # exact frame index
                    Mcross.tStart = t  # local t and not account for scr refresh
                    Mcross.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(Mcross, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'Mcross.started')
                    # update status
                    Mcross.status = STARTED
                    Mcross.setAutoDraw(True)
                
                # if Mcross is active this frame...
                if Mcross.status == STARTED:
                    # update params
                    pass
                
                # if Mcross is stopping this frame...
                if Mcross.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > Mcross.tStartRefresh + 1-frameTolerance:
                        # keep track of stop time/frame for later
                        Mcross.tStop = t  # not accounting for scr refresh
                        Mcross.tStopRefresh = tThisFlipGlobal  # on global time
                        Mcross.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'Mcross.stopped')
                        # update status
                        Mcross.status = FINISHED
                        Mcross.setAutoDraw(False)
                
                # check for quit (typically the Esc key)
                if defaultKeyboard.getKeys(keyList=["escape"]):
                    thisExp.status = FINISHED
                if thisExp.status == FINISHED or endExpNow:
                    endExperiment(thisExp, win=win)
                    return
                # pause experiment here if requested
                if thisExp.status == PAUSED:
                    pauseExperiment(
                        thisExp=thisExp, 
                        win=win, 
                        timers=[routineTimer], 
                        playbackComponents=[]
                    )
                    # skip the frame we paused on
                    continue
                
                # check if all components have finished
                if not continueRoutine:  # a component has requested a forced-end of Routine
                    feedback.forceEnded = routineForceEnded = True
                    break
                continueRoutine = False  # will revert to True if at least one component still running
                for thisComponent in feedback.components:
                    if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                        continueRoutine = True
                        break  # at least one component has not yet finished
                
                # refresh the screen
                if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                    win.flip()
            
            # --- Ending Routine "feedback" ---
            for thisComponent in feedback.components:
                if hasattr(thisComponent, "setAutoDraw"):
                    thisComponent.setAutoDraw(False)
            # store stop times for feedback
            feedback.tStop = globalClock.getTime(format='float')
            feedback.tStopRefresh = tThisFlipGlobal
            thisExp.addData('feedback.stopped', feedback.tStop)
            # Run 'End Routine' code from fb_code
            thisExp.addData('correct', correct)
            thisExp.addData('outcome', outcome)
            thisExp.addData('trial_key', trial_key)
            # the Routine "feedback" was not non-slip safe, so reset the non-slip timer
            routineTimer.reset()
            
            # --- Prepare to start Routine "ITI" ---
            # create an object to store info about Routine ITI
            ITI = data.Routine(
                name='ITI',
                components=[],
            )
            ITI.status = NOT_STARTED
            continueRoutine = True
            # update component parameters for each repeat
            # store start times for ITI
            ITI.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
            ITI.tStart = globalClock.getTime(format='float')
            ITI.status = STARTED
            thisExp.addData('ITI.started', ITI.tStart)
            ITI.maxDuration = 1
            # keep track of which components have finished
            ITIComponents = ITI.components
            for thisComponent in ITI.components:
                thisComponent.tStart = None
                thisComponent.tStop = None
                thisComponent.tStartRefresh = None
                thisComponent.tStopRefresh = None
                if hasattr(thisComponent, 'status'):
                    thisComponent.status = NOT_STARTED
            # reset timers
            t = 0
            _timeToFirstFrame = win.getFutureFlipTime(clock="now")
            frameN = -1
            
            # --- Run Routine "ITI" ---
            # if trial has changed, end Routine now
            if isinstance(trials, data.TrialHandler2) and thisTrial.thisN != trials.thisTrial.thisN:
                continueRoutine = False
            ITI.forceEnded = routineForceEnded = not continueRoutine
            while continueRoutine and routineTimer.getTime() < 1.0:
                # get current time
                t = routineTimer.getTime()
                tThisFlip = win.getFutureFlipTime(clock=routineTimer)
                tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                # update/draw components on each frame
                # is it time to end the Routine? (based on local clock)
                if tThisFlip > ITI.maxDuration-frameTolerance:
                    ITI.maxDurationReached = True
                    continueRoutine = False
                
                # check for quit (typically the Esc key)
                if defaultKeyboard.getKeys(keyList=["escape"]):
                    thisExp.status = FINISHED
                if thisExp.status == FINISHED or endExpNow:
                    endExperiment(thisExp, win=win)
                    return
                # pause experiment here if requested
                if thisExp.status == PAUSED:
                    pauseExperiment(
                        thisExp=thisExp, 
                        win=win, 
                        timers=[routineTimer], 
                        playbackComponents=[]
                    )
                    # skip the frame we paused on
                    continue
                
                # check if all components have finished
                if not continueRoutine:  # a component has requested a forced-end of Routine
                    ITI.forceEnded = routineForceEnded = True
                    break
                continueRoutine = False  # will revert to True if at least one component still running
                for thisComponent in ITI.components:
                    if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                        continueRoutine = True
                        break  # at least one component has not yet finished
                
                # refresh the screen
                if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                    win.flip()
            
            # --- Ending Routine "ITI" ---
            for thisComponent in ITI.components:
                if hasattr(thisComponent, "setAutoDraw"):
                    thisComponent.setAutoDraw(False)
            # store stop times for ITI
            ITI.tStop = globalClock.getTime(format='float')
            ITI.tStopRefresh = tThisFlipGlobal
            thisExp.addData('ITI.stopped', ITI.tStop)
            # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
            if ITI.maxDurationReached:
                routineTimer.addTime(-ITI.maxDuration)
            elif ITI.forceEnded:
                routineTimer.reset()
            else:
                routineTimer.addTime(-1.000000)
            thisExp.nextEntry()
            
        # completed 1.0 repeats of 'trials'
        
        if thisSession is not None:
            # if running in a Session with a Liaison client, send data up to now
            thisSession.sendExperimentData()
        
        # --- Prepare to start Routine "end_block" ---
        # create an object to store info about Routine end_block
        end_block = data.Routine(
            name='end_block',
            components=[block_end_text],
        )
        end_block.status = NOT_STARTED
        continueRoutine = True
        # update component parameters for each repeat
        # Run 'Begin Routine' code from block_end_code
        block_end_txt = f'You gained {block_outcome} big coins,\nand a bonus of {int(np.ceil(block_bonus))} small coins'
        block_end_text.setText(block_end_txt)
        # store start times for end_block
        end_block.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
        end_block.tStart = globalClock.getTime(format='float')
        end_block.status = STARTED
        thisExp.addData('end_block.started', end_block.tStart)
        end_block.maxDuration = None
        # keep track of which components have finished
        end_blockComponents = end_block.components
        for thisComponent in end_block.components:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "end_block" ---
        # if trial has changed, end Routine now
        if isinstance(blocks, data.TrialHandler2) and thisBlock.thisN != blocks.thisTrial.thisN:
            continueRoutine = False
        end_block.forceEnded = routineForceEnded = not continueRoutine
        while continueRoutine and routineTimer.getTime() < 4.0:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *block_end_text* updates
            
            # if block_end_text is starting this frame...
            if block_end_text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                block_end_text.frameNStart = frameN  # exact frame index
                block_end_text.tStart = t  # local t and not account for scr refresh
                block_end_text.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(block_end_text, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'block_end_text.started')
                # update status
                block_end_text.status = STARTED
                block_end_text.setAutoDraw(True)
            
            # if block_end_text is active this frame...
            if block_end_text.status == STARTED:
                # update params
                pass
            
            # if block_end_text is stopping this frame...
            if block_end_text.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > block_end_text.tStartRefresh + 4-frameTolerance:
                    # keep track of stop time/frame for later
                    block_end_text.tStop = t  # not accounting for scr refresh
                    block_end_text.tStopRefresh = tThisFlipGlobal  # on global time
                    block_end_text.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'block_end_text.stopped')
                    # update status
                    block_end_text.status = FINISHED
                    block_end_text.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, win=win)
                return
            # pause experiment here if requested
            if thisExp.status == PAUSED:
                pauseExperiment(
                    thisExp=thisExp, 
                    win=win, 
                    timers=[routineTimer], 
                    playbackComponents=[]
                )
                # skip the frame we paused on
                continue
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                end_block.forceEnded = routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in end_block.components:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "end_block" ---
        for thisComponent in end_block.components:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # store stop times for end_block
        end_block.tStop = globalClock.getTime(format='float')
        end_block.tStopRefresh = tThisFlipGlobal
        thisExp.addData('end_block.stopped', end_block.tStop)
        # Run 'End Routine' code from block_end_code
        thisExp.addData('block_outcome', block_outcome)
        thisExp.addData('block_bonus', block_bonus)
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if end_block.maxDurationReached:
            routineTimer.addTime(-end_block.maxDuration)
        elif end_block.forceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-4.000000)
    # completed 1.0 repeats of 'blocks'
    
    
    # mark experiment as finished
    endExperiment(thisExp, win=win)


def saveData(thisExp):
    """
    Save data from this experiment
    
    Parameters
    ==========
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    """
    filename = thisExp.dataFileName
    # these shouldn't be strictly necessary (should auto-save)
    thisExp.saveAsWideText(filename + '.csv', delim='auto')
    thisExp.saveAsPickle(filename)


def endExperiment(thisExp, win=None):
    """
    End this experiment, performing final shut down operations.
    
    This function does NOT close the window or end the Python process - use `quit` for this.
    
    Parameters
    ==========
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    win : psychopy.visual.Window
        Window for this experiment.
    """
    if win is not None:
        # remove autodraw from all current components
        win.clearAutoDraw()
        # Flip one final time so any remaining win.callOnFlip() 
        # and win.timeOnFlip() tasks get executed
        win.flip()
    # return console logger level to WARNING
    logging.console.setLevel(logging.WARNING)
    # mark experiment handler as finished
    thisExp.status = FINISHED


def quit(thisExp, win=None, thisSession=None):
    """
    Fully quit, closing the window and ending the Python process.
    
    Parameters
    ==========
    win : psychopy.visual.Window
        Window to close.
    thisSession : psychopy.session.Session or None
        Handle of the Session object this experiment is being run from, if any.
    """
    thisExp.abort()  # or data files will save again on exit
    # make sure everything is closed down
    if win is not None:
        # Flip one final time so any remaining win.callOnFlip() 
        # and win.timeOnFlip() tasks get executed before quitting
        win.flip()
        win.close()
    if thisSession is not None:
        thisSession.stop()
    # terminate Python process
    core.quit()


# if running this experiment as a script...
if __name__ == '__main__':
    # call all functions in order
    expInfo = showExpInfoDlg(expInfo=expInfo)
    thisExp = setupData(expInfo=expInfo)
    logFile = setupLogging(filename=thisExp.dataFileName)
    win = setupWindow(expInfo=expInfo)
    setupDevices(expInfo=expInfo, thisExp=thisExp, win=win)
    run(
        expInfo=expInfo, 
        thisExp=thisExp, 
        win=win,
        globalClock='float'
    )
    saveData(thisExp=thisExp)
    quit(thisExp=thisExp, win=win)
