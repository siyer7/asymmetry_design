#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v2024.2.4),
    on Sat Mar 15 20:53:47 2025
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
expName = 'asymmetry_practice'  # from the Builder filename that created this script
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
    filename = u'results/%s_subj-%s_order%s_%s' % (expName, expInfo['subj'], expInfo['sess_type'], expInfo['date'])
    # make sure filename is relative to dataDir
    if os.path.isabs(filename):
        dataDir = os.path.commonprefix([dataDir, filename])
        filename = os.path.relpath(filename, dataDir)
    
    # an ExperimentHandler isn't essential but helps with data saving
    thisExp = data.ExperimentHandler(
        name=expName, version='',
        extraInfo=expInfo, runtimeInfo=None,
        originPath='/Users/f0064z8/Library/CloudStorage/GoogleDrive-si2442@columbia.edu/My Drive/research/asymmetry_design/asymmetry_tutorial_lastrun.py',
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
            size=_winSize, fullscr=_fullScr, screen=1,
            winType='pyglet', allowGUI=True, allowStencil=False,
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
        expInfo['frameRate'] = 120
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
    if deviceManager.getDevice('intro_resp') is None:
        # initialise intro_resp
        intro_resp = deviceManager.addDevice(
            deviceClass='keyboard',
            deviceName='intro_resp',
        )
    if deviceManager.getDevice('expl_resp1') is None:
        # initialise expl_resp1
        expl_resp1 = deviceManager.addDevice(
            deviceClass='keyboard',
            deviceName='expl_resp1',
        )
    if deviceManager.getDevice('expl_resp') is None:
        # initialise expl_resp
        expl_resp = deviceManager.addDevice(
            deviceClass='keyboard',
            deviceName='expl_resp',
        )
    if deviceManager.getDevice('base_resp') is None:
        # initialise base_resp
        base_resp = deviceManager.addDevice(
            deviceClass='keyboard',
            deviceName='base_resp',
        )
    if deviceManager.getDevice('stim_resp') is None:
        # initialise stim_resp
        stim_resp = deviceManager.addDevice(
            deviceClass='keyboard',
            deviceName='stim_resp',
        )
    if deviceManager.getDevice('delay_resp') is None:
        # initialise delay_resp
        delay_resp = deviceManager.addDevice(
            deviceClass='keyboard',
            deviceName='delay_resp',
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
    if deviceManager.getDevice('feedback_resp') is None:
        # initialise feedback_resp
        feedback_resp = deviceManager.addDevice(
            deviceClass='keyboard',
            deviceName='feedback_resp',
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
    
    # --- Initialize components for Routine "intro" ---
    # Run 'Begin Experiment' code from reward_struct_code
    if expInfo['sess_type'] in ['A','C']:
        img_good_file = 'input_data/stims/Shape_4.00_2.00.png'
        img_bad_file = 'input_data/stims/Shape_0.00_2.00.png'
    
    else:# expInfo['sess_type'] in ['B','D']:
        img_good_file = 'input_data/stims/Shape_0.00_2.00.png'
        img_bad_file = 'input_data/stims/Shape_4.00_2.00.png'
        
        
    correct_pos, wrong_pos = .2, .4
    intro_text = visual.TextStim(win=win, name='intro_text',
        text='Shape identification task\n\nMaximize reward and bonus',
        font='Arial',
        pos=(0, 0), draggable=False, height=0.05, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-1.0);
    intro_resp = keyboard.Keyboard(deviceName='intro_resp')
    
    # --- Initialize components for Routine "expl1" ---
    expl_resp1 = keyboard.Keyboard(deviceName='expl_resp1')
    good1 = visual.TextStim(win=win, name='good1',
        text='',
        font='Arial',
        pos=[-.4,.2], draggable=False, height=0.05, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-1.0);
    img_good1 = visual.ImageStim(
        win=win,
        name='img_good1', 
        image=img_good_file, mask=None, anchor='center',
        ori=0.0, pos=[0,.2], draggable=False, size=[.4,.4],
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-2.0)
    correct1 = visual.TextStim(win=win, name='correct1',
        text='',
        font='Arial',
        pos=[correct_pos,.3], draggable=False, height=0.05, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-3.0);
    wrong1 = visual.TextStim(win=win, name='wrong1',
        text='',
        font='Arial',
        pos=[wrong_pos,.3], draggable=False, height=0.05, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-4.0);
    highR1 = visual.TextStim(win=win, name='highR1',
        text='',
        font='Arial',
        pos=[correct_pos,.2], draggable=False, height=0.05, wrapWidth=None, ori=0.0, 
        color=[-1.0000, 0.5000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-5.0);
    lowP1 = visual.TextStim(win=win, name='lowP1',
        text='',
        font='Arial',
        pos=[wrong_pos,.2], draggable=False, height=0.05, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-6.0);
    
    # --- Initialize components for Routine "expl2" ---
    expl_resp = keyboard.Keyboard(deviceName='expl_resp')
    good2 = visual.TextStim(win=win, name='good2',
        text='',
        font='Arial',
        pos=[-.4,.2], draggable=False, height=0.05, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-1.0);
    img_good2 = visual.ImageStim(
        win=win,
        name='img_good2', 
        image=img_good_file, mask=None, anchor='center',
        ori=0.0, pos=[0,.2], draggable=False, size=[.4,.4],
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-2.0)
    correct2 = visual.TextStim(win=win, name='correct2',
        text='',
        font='Arial',
        pos=[correct_pos,.3], draggable=False, height=0.05, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-3.0);
    highR2 = visual.TextStim(win=win, name='highR2',
        text='',
        font='Arial',
        pos=[correct_pos,.2], draggable=False, height=0.05, wrapWidth=None, ori=0.0, 
        color=[-1.0000, 0.5000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-4.0);
    lowP2 = visual.TextStim(win=win, name='lowP2',
        text='',
        font='Arial',
        pos=[wrong_pos,.2], draggable=False, height=0.05, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-5.0);
    bad = visual.TextStim(win=win, name='bad',
        text='',
        font='Arial',
        pos=[-.4,-.2], draggable=False, height=0.05, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-6.0);
    img_bad = visual.ImageStim(
        win=win,
        name='img_bad', 
        image=img_bad_file, mask=None, anchor='center',
        ori=0.0, pos=[0,-.2], draggable=False, size=[.4,.4],
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-7.0)
    wrong = visual.TextStim(win=win, name='wrong',
        text='',
        font='Arial',
        pos=[wrong_pos,.3], draggable=False, height=0.05, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-8.0);
    lowR = visual.TextStim(win=win, name='lowR',
        text='',
        font='Arial',
        pos=[correct_pos,-.2], draggable=False, height=0.05, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-9.0);
    highP = visual.TextStim(win=win, name='highP',
        text='',
        font='Arial',
        pos=[wrong_pos,-.2], draggable=False, height=0.05, wrapWidth=None, ori=0.0, 
        color=[1.0000, 0.5000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-10.0);
    
    # --- Initialize components for Routine "baseline" ---
    ISI1 = clock.StaticPeriod(win=win, screenHz=expInfo['frameRate'], name='ISI1')
    base_resp = keyboard.Keyboard(deviceName='base_resp')
    
    # --- Initialize components for Routine "stim" ---
    target_stim = visual.ImageStim(
        win=win,
        name='target_stim', 
        image='default.png', mask=None, anchor='center',
        ori=0.0, pos=(0, 0), draggable=False, size=[.6,.6],
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-1.0)
    stim_resp = keyboard.Keyboard(deviceName='stim_resp')
    
    # --- Initialize components for Routine "delay" ---
    ISI2 = clock.StaticPeriod(win=win, screenHz=expInfo['frameRate'], name='ISI2')
    delay_resp = keyboard.Keyboard(deviceName='delay_resp')
    
    # --- Initialize components for Routine "task" ---
    # Run 'Begin Experiment' code from rand_divider_slider
    leftPressed, rightPressed, marker_moved = .0, .0, .0
    positions = []
    marker_move = .004
    debug_task_txt = ''
    img1 = visual.ImageStim(
        win=win,
        name='img1', 
        image='default.png', mask=None, anchor='center',
        ori=0.0, pos=[-.6, 0], draggable=False, size=[.45,.45],
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-1.0)
    img2 = visual.ImageStim(
        win=win,
        name='img2', 
        image='default.png', mask=None, anchor='center',
        ori=0.0, pos=[.6, 0], draggable=False, size=[.45,.45],
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-2.0)
    divider_line = visual.Line(
        win=win, name='divider_line',
        size=[.4,0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=5.0,
        colorSpace='rgb', lineColor=[1.0000, 1.0000, 1.0000], fillColor=[1.0000, 1.0000, 1.0000],
        opacity=None, depth=-3.0, interpolate=True)
    slider_line = visual.Line(
        win=win, name='slider_line',
        size=[.8,.0],
        ori=0.0, pos=[0, 0], draggable=False, anchor='center',
        lineWidth=5.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=1.0, depth=-4.0, interpolate=True)
    marker = visual.Slider(win=win, name='marker',
        startValue=None, size=[.8,.05], pos=[0,0], units=win.units,
        labels=None, ticks=(-.4,-.3,-.2,-.1,.0,.1,.2,.3,.4), granularity=0.001,
        style='slider', styleTweaks=('labels45',), opacity=1.0,
        labelColor=None, markerColor=[-1.0000, -1.0000, -1.0000], lineColor=None, colorSpace='rgb',
        font='Open Sans', labelHeight=0.05,
        flip=False, ori=0.0, depth=-5, readOnly=False)
    slider_resp = keyboard.Keyboard(deviceName='slider_resp')
    task_text = visual.TextStim(win=win, name='task_text',
        text='',
        font='Arial',
        pos=[0,0], draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-7.0);
    submit_resp = keyboard.Keyboard(deviceName='submit_resp')
    
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
    feedback_resp = keyboard.Keyboard(deviceName='feedback_resp')
    
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
    
    # --- Prepare to start Routine "intro" ---
    # create an object to store info about Routine intro
    intro = data.Routine(
        name='intro',
        components=[intro_text, intro_resp],
    )
    intro.status = NOT_STARTED
    continueRoutine = True
    # update component parameters for each repeat
    # create starting attributes for intro_resp
    intro_resp.keys = []
    intro_resp.rt = []
    _intro_resp_allKeys = []
    # store start times for intro
    intro.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
    intro.tStart = globalClock.getTime(format='float')
    intro.status = STARTED
    thisExp.addData('intro.started', intro.tStart)
    intro.maxDuration = None
    # keep track of which components have finished
    introComponents = intro.components
    for thisComponent in intro.components:
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
    
    # --- Run Routine "intro" ---
    intro.forceEnded = routineForceEnded = not continueRoutine
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *intro_text* updates
        
        # if intro_text is starting this frame...
        if intro_text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            intro_text.frameNStart = frameN  # exact frame index
            intro_text.tStart = t  # local t and not account for scr refresh
            intro_text.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(intro_text, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'intro_text.started')
            # update status
            intro_text.status = STARTED
            intro_text.setAutoDraw(True)
        
        # if intro_text is active this frame...
        if intro_text.status == STARTED:
            # update params
            pass
        
        # *intro_resp* updates
        waitOnFlip = False
        
        # if intro_resp is starting this frame...
        if intro_resp.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            intro_resp.frameNStart = frameN  # exact frame index
            intro_resp.tStart = t  # local t and not account for scr refresh
            intro_resp.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(intro_resp, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'intro_resp.started')
            # update status
            intro_resp.status = STARTED
            # keyboard checking is just starting
            waitOnFlip = True
            win.callOnFlip(intro_resp.clock.reset)  # t=0 on next screen flip
            win.callOnFlip(intro_resp.clearEvents, eventType='keyboard')  # clear events on next screen flip
        if intro_resp.status == STARTED and not waitOnFlip:
            theseKeys = intro_resp.getKeys(keyList=['up'], ignoreKeys=["escape"], waitRelease=False)
            _intro_resp_allKeys.extend(theseKeys)
            if len(_intro_resp_allKeys):
                intro_resp.keys = _intro_resp_allKeys[-1].name  # just the last key pressed
                intro_resp.rt = _intro_resp_allKeys[-1].rt
                intro_resp.duration = _intro_resp_allKeys[-1].duration
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
            intro.forceEnded = routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in intro.components:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "intro" ---
    for thisComponent in intro.components:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # store stop times for intro
    intro.tStop = globalClock.getTime(format='float')
    intro.tStopRefresh = tThisFlipGlobal
    thisExp.addData('intro.stopped', intro.tStop)
    # check responses
    if intro_resp.keys in ['', [], None]:  # No response was made
        intro_resp.keys = None
    thisExp.addData('intro_resp.keys',intro_resp.keys)
    if intro_resp.keys != None:  # we had a response
        thisExp.addData('intro_resp.rt', intro_resp.rt)
        thisExp.addData('intro_resp.duration', intro_resp.duration)
    thisExp.nextEntry()
    # the Routine "intro" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # --- Prepare to start Routine "expl1" ---
    # create an object to store info about Routine expl1
    expl1 = data.Routine(
        name='expl1',
        components=[expl_resp1, good1, img_good1, correct1, wrong1, highR1, lowP1],
    )
    expl1.status = NOT_STARTED
    continueRoutine = True
    # update component parameters for each repeat
    # create starting attributes for expl_resp1
    expl_resp1.keys = []
    expl_resp1.rt = []
    _expl_resp1_allKeys = []
    good1.setText('good')
    correct1.setText('correct')
    wrong1.setText('wrong')
    highR1.setText('2')
    lowP1.setText('-1')
    # store start times for expl1
    expl1.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
    expl1.tStart = globalClock.getTime(format='float')
    expl1.status = STARTED
    thisExp.addData('expl1.started', expl1.tStart)
    expl1.maxDuration = None
    # keep track of which components have finished
    expl1Components = expl1.components
    for thisComponent in expl1.components:
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
    
    # --- Run Routine "expl1" ---
    expl1.forceEnded = routineForceEnded = not continueRoutine
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *expl_resp1* updates
        waitOnFlip = False
        
        # if expl_resp1 is starting this frame...
        if expl_resp1.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            expl_resp1.frameNStart = frameN  # exact frame index
            expl_resp1.tStart = t  # local t and not account for scr refresh
            expl_resp1.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(expl_resp1, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'expl_resp1.started')
            # update status
            expl_resp1.status = STARTED
            # keyboard checking is just starting
            waitOnFlip = True
            win.callOnFlip(expl_resp1.clock.reset)  # t=0 on next screen flip
            win.callOnFlip(expl_resp1.clearEvents, eventType='keyboard')  # clear events on next screen flip
        if expl_resp1.status == STARTED and not waitOnFlip:
            theseKeys = expl_resp1.getKeys(keyList=['up'], ignoreKeys=["escape"], waitRelease=False)
            _expl_resp1_allKeys.extend(theseKeys)
            if len(_expl_resp1_allKeys):
                expl_resp1.keys = _expl_resp1_allKeys[-1].name  # just the last key pressed
                expl_resp1.rt = _expl_resp1_allKeys[-1].rt
                expl_resp1.duration = _expl_resp1_allKeys[-1].duration
                # a response ends the routine
                continueRoutine = False
        
        # *good1* updates
        
        # if good1 is starting this frame...
        if good1.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            good1.frameNStart = frameN  # exact frame index
            good1.tStart = t  # local t and not account for scr refresh
            good1.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(good1, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'good1.started')
            # update status
            good1.status = STARTED
            good1.setAutoDraw(True)
        
        # if good1 is active this frame...
        if good1.status == STARTED:
            # update params
            pass
        
        # *img_good1* updates
        
        # if img_good1 is starting this frame...
        if img_good1.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            img_good1.frameNStart = frameN  # exact frame index
            img_good1.tStart = t  # local t and not account for scr refresh
            img_good1.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(img_good1, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'img_good1.started')
            # update status
            img_good1.status = STARTED
            img_good1.setAutoDraw(True)
        
        # if img_good1 is active this frame...
        if img_good1.status == STARTED:
            # update params
            pass
        
        # *correct1* updates
        
        # if correct1 is starting this frame...
        if correct1.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            correct1.frameNStart = frameN  # exact frame index
            correct1.tStart = t  # local t and not account for scr refresh
            correct1.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(correct1, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'correct1.started')
            # update status
            correct1.status = STARTED
            correct1.setAutoDraw(True)
        
        # if correct1 is active this frame...
        if correct1.status == STARTED:
            # update params
            pass
        
        # *wrong1* updates
        
        # if wrong1 is starting this frame...
        if wrong1.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            wrong1.frameNStart = frameN  # exact frame index
            wrong1.tStart = t  # local t and not account for scr refresh
            wrong1.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(wrong1, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'wrong1.started')
            # update status
            wrong1.status = STARTED
            wrong1.setAutoDraw(True)
        
        # if wrong1 is active this frame...
        if wrong1.status == STARTED:
            # update params
            pass
        
        # *highR1* updates
        
        # if highR1 is starting this frame...
        if highR1.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            highR1.frameNStart = frameN  # exact frame index
            highR1.tStart = t  # local t and not account for scr refresh
            highR1.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(highR1, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'highR1.started')
            # update status
            highR1.status = STARTED
            highR1.setAutoDraw(True)
        
        # if highR1 is active this frame...
        if highR1.status == STARTED:
            # update params
            pass
        
        # *lowP1* updates
        
        # if lowP1 is starting this frame...
        if lowP1.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            lowP1.frameNStart = frameN  # exact frame index
            lowP1.tStart = t  # local t and not account for scr refresh
            lowP1.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(lowP1, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'lowP1.started')
            # update status
            lowP1.status = STARTED
            lowP1.setAutoDraw(True)
        
        # if lowP1 is active this frame...
        if lowP1.status == STARTED:
            # update params
            pass
        
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
            expl1.forceEnded = routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in expl1.components:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "expl1" ---
    for thisComponent in expl1.components:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # store stop times for expl1
    expl1.tStop = globalClock.getTime(format='float')
    expl1.tStopRefresh = tThisFlipGlobal
    thisExp.addData('expl1.stopped', expl1.tStop)
    # check responses
    if expl_resp1.keys in ['', [], None]:  # No response was made
        expl_resp1.keys = None
    thisExp.addData('expl_resp1.keys',expl_resp1.keys)
    if expl_resp1.keys != None:  # we had a response
        thisExp.addData('expl_resp1.rt', expl_resp1.rt)
        thisExp.addData('expl_resp1.duration', expl_resp1.duration)
    thisExp.nextEntry()
    # the Routine "expl1" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # --- Prepare to start Routine "expl2" ---
    # create an object to store info about Routine expl2
    expl2 = data.Routine(
        name='expl2',
        components=[expl_resp, good2, img_good2, correct2, highR2, lowP2, bad, img_bad, wrong, lowR, highP],
    )
    expl2.status = NOT_STARTED
    continueRoutine = True
    # update component parameters for each repeat
    # create starting attributes for expl_resp
    expl_resp.keys = []
    expl_resp.rt = []
    _expl_resp_allKeys = []
    good2.setText('good')
    correct2.setText('correct')
    highR2.setText('2')
    lowP2.setText('-1')
    bad.setText('bad')
    wrong.setText('wrong')
    lowR.setText('1')
    highP.setText('-2')
    # store start times for expl2
    expl2.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
    expl2.tStart = globalClock.getTime(format='float')
    expl2.status = STARTED
    thisExp.addData('expl2.started', expl2.tStart)
    expl2.maxDuration = None
    # keep track of which components have finished
    expl2Components = expl2.components
    for thisComponent in expl2.components:
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
    
    # --- Run Routine "expl2" ---
    expl2.forceEnded = routineForceEnded = not continueRoutine
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *expl_resp* updates
        waitOnFlip = False
        
        # if expl_resp is starting this frame...
        if expl_resp.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            expl_resp.frameNStart = frameN  # exact frame index
            expl_resp.tStart = t  # local t and not account for scr refresh
            expl_resp.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(expl_resp, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'expl_resp.started')
            # update status
            expl_resp.status = STARTED
            # keyboard checking is just starting
            waitOnFlip = True
            win.callOnFlip(expl_resp.clock.reset)  # t=0 on next screen flip
            win.callOnFlip(expl_resp.clearEvents, eventType='keyboard')  # clear events on next screen flip
        if expl_resp.status == STARTED and not waitOnFlip:
            theseKeys = expl_resp.getKeys(keyList=['up'], ignoreKeys=["escape"], waitRelease=False)
            _expl_resp_allKeys.extend(theseKeys)
            if len(_expl_resp_allKeys):
                expl_resp.keys = _expl_resp_allKeys[-1].name  # just the last key pressed
                expl_resp.rt = _expl_resp_allKeys[-1].rt
                expl_resp.duration = _expl_resp_allKeys[-1].duration
                # a response ends the routine
                continueRoutine = False
        
        # *good2* updates
        
        # if good2 is starting this frame...
        if good2.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            good2.frameNStart = frameN  # exact frame index
            good2.tStart = t  # local t and not account for scr refresh
            good2.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(good2, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'good2.started')
            # update status
            good2.status = STARTED
            good2.setAutoDraw(True)
        
        # if good2 is active this frame...
        if good2.status == STARTED:
            # update params
            pass
        
        # *img_good2* updates
        
        # if img_good2 is starting this frame...
        if img_good2.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            img_good2.frameNStart = frameN  # exact frame index
            img_good2.tStart = t  # local t and not account for scr refresh
            img_good2.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(img_good2, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'img_good2.started')
            # update status
            img_good2.status = STARTED
            img_good2.setAutoDraw(True)
        
        # if img_good2 is active this frame...
        if img_good2.status == STARTED:
            # update params
            pass
        
        # *correct2* updates
        
        # if correct2 is starting this frame...
        if correct2.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            correct2.frameNStart = frameN  # exact frame index
            correct2.tStart = t  # local t and not account for scr refresh
            correct2.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(correct2, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'correct2.started')
            # update status
            correct2.status = STARTED
            correct2.setAutoDraw(True)
        
        # if correct2 is active this frame...
        if correct2.status == STARTED:
            # update params
            pass
        
        # *highR2* updates
        
        # if highR2 is starting this frame...
        if highR2.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            highR2.frameNStart = frameN  # exact frame index
            highR2.tStart = t  # local t and not account for scr refresh
            highR2.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(highR2, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'highR2.started')
            # update status
            highR2.status = STARTED
            highR2.setAutoDraw(True)
        
        # if highR2 is active this frame...
        if highR2.status == STARTED:
            # update params
            pass
        
        # *lowP2* updates
        
        # if lowP2 is starting this frame...
        if lowP2.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            lowP2.frameNStart = frameN  # exact frame index
            lowP2.tStart = t  # local t and not account for scr refresh
            lowP2.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(lowP2, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'lowP2.started')
            # update status
            lowP2.status = STARTED
            lowP2.setAutoDraw(True)
        
        # if lowP2 is active this frame...
        if lowP2.status == STARTED:
            # update params
            pass
        
        # *bad* updates
        
        # if bad is starting this frame...
        if bad.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            bad.frameNStart = frameN  # exact frame index
            bad.tStart = t  # local t and not account for scr refresh
            bad.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(bad, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'bad.started')
            # update status
            bad.status = STARTED
            bad.setAutoDraw(True)
        
        # if bad is active this frame...
        if bad.status == STARTED:
            # update params
            pass
        
        # *img_bad* updates
        
        # if img_bad is starting this frame...
        if img_bad.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            img_bad.frameNStart = frameN  # exact frame index
            img_bad.tStart = t  # local t and not account for scr refresh
            img_bad.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(img_bad, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'img_bad.started')
            # update status
            img_bad.status = STARTED
            img_bad.setAutoDraw(True)
        
        # if img_bad is active this frame...
        if img_bad.status == STARTED:
            # update params
            pass
        
        # *wrong* updates
        
        # if wrong is starting this frame...
        if wrong.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            wrong.frameNStart = frameN  # exact frame index
            wrong.tStart = t  # local t and not account for scr refresh
            wrong.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(wrong, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'wrong.started')
            # update status
            wrong.status = STARTED
            wrong.setAutoDraw(True)
        
        # if wrong is active this frame...
        if wrong.status == STARTED:
            # update params
            pass
        
        # *lowR* updates
        
        # if lowR is starting this frame...
        if lowR.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            lowR.frameNStart = frameN  # exact frame index
            lowR.tStart = t  # local t and not account for scr refresh
            lowR.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(lowR, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'lowR.started')
            # update status
            lowR.status = STARTED
            lowR.setAutoDraw(True)
        
        # if lowR is active this frame...
        if lowR.status == STARTED:
            # update params
            pass
        
        # *highP* updates
        
        # if highP is starting this frame...
        if highP.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            highP.frameNStart = frameN  # exact frame index
            highP.tStart = t  # local t and not account for scr refresh
            highP.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(highP, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'highP.started')
            # update status
            highP.status = STARTED
            highP.setAutoDraw(True)
        
        # if highP is active this frame...
        if highP.status == STARTED:
            # update params
            pass
        
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
            expl2.forceEnded = routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in expl2.components:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "expl2" ---
    for thisComponent in expl2.components:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # store stop times for expl2
    expl2.tStop = globalClock.getTime(format='float')
    expl2.tStopRefresh = tThisFlipGlobal
    thisExp.addData('expl2.stopped', expl2.tStop)
    # check responses
    if expl_resp.keys in ['', [], None]:  # No response was made
        expl_resp.keys = None
    thisExp.addData('expl_resp.keys',expl_resp.keys)
    if expl_resp.keys != None:  # we had a response
        thisExp.addData('expl_resp.rt', expl_resp.rt)
        thisExp.addData('expl_resp.duration', expl_resp.duration)
    thisExp.nextEntry()
    # the Routine "expl2" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # set up handler to look after randomisation of conditions etc
    trials = data.TrialHandler2(
        name='trials',
        nReps=1.0, 
        method='sequential', 
        extraInfo=expInfo, 
        originPath=-1, 
        trialList=data.importConditions('input_data/tutorial.csv'), 
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
        
        # --- Prepare to start Routine "baseline" ---
        # create an object to store info about Routine baseline
        baseline = data.Routine(
            name='baseline',
            components=[ISI1, base_resp],
        )
        baseline.status = NOT_STARTED
        continueRoutine = True
        # update component parameters for each repeat
        # create starting attributes for base_resp
        base_resp.keys = []
        base_resp.rt = []
        _base_resp_allKeys = []
        # store start times for baseline
        baseline.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
        baseline.tStart = globalClock.getTime(format='float')
        baseline.status = STARTED
        thisExp.addData('baseline.started', baseline.tStart)
        # keep track of which components have finished
        baselineComponents = baseline.components
        for thisComponent in baseline.components:
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
        
        # --- Run Routine "baseline" ---
        # if trial has changed, end Routine now
        if isinstance(trials, data.TrialHandler2) and thisTrial.thisN != trials.thisTrial.thisN:
            continueRoutine = False
        baseline.forceEnded = routineForceEnded = not continueRoutine
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *base_resp* updates
            waitOnFlip = False
            
            # if base_resp is starting this frame...
            if base_resp.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                base_resp.frameNStart = frameN  # exact frame index
                base_resp.tStart = t  # local t and not account for scr refresh
                base_resp.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(base_resp, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'base_resp.started')
                # update status
                base_resp.status = STARTED
                # keyboard checking is just starting
                waitOnFlip = True
                win.callOnFlip(base_resp.clock.reset)  # t=0 on next screen flip
                win.callOnFlip(base_resp.clearEvents, eventType='keyboard')  # clear events on next screen flip
            if base_resp.status == STARTED and not waitOnFlip:
                theseKeys = base_resp.getKeys(keyList=['up'], ignoreKeys=["escape"], waitRelease=False)
                _base_resp_allKeys.extend(theseKeys)
                if len(_base_resp_allKeys):
                    base_resp.keys = _base_resp_allKeys[-1].name  # just the last key pressed
                    base_resp.rt = _base_resp_allKeys[-1].rt
                    base_resp.duration = _base_resp_allKeys[-1].duration
                    # a response ends the routine
                    continueRoutine = False
            # *ISI1* period
            
            # if ISI1 is starting this frame...
            if ISI1.status == NOT_STARTED and frameN >= 0:
                # keep track of start time/frame for later
                ISI1.frameNStart = frameN  # exact frame index
                ISI1.tStart = t  # local t and not account for scr refresh
                ISI1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(ISI1, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.addData('ISI1.started', t)
                # update status
                ISI1.status = STARTED
                ISI1.start(1*frameDur)
            elif ISI1.status == STARTED:  # one frame should pass before updating params and completing
                # Updating other components during *ISI1*
                target_stim.setImage(target_file)
                # Component updates done
                ISI1.complete()  # finish the static period
                ISI1.tStop = ISI1.tStart + 1*frameDur  # record stop time
            
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
                baseline.forceEnded = routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in baseline.components:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "baseline" ---
        for thisComponent in baseline.components:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # store stop times for baseline
        baseline.tStop = globalClock.getTime(format='float')
        baseline.tStopRefresh = tThisFlipGlobal
        thisExp.addData('baseline.stopped', baseline.tStop)
        # check responses
        if base_resp.keys in ['', [], None]:  # No response was made
            base_resp.keys = None
        trials.addData('base_resp.keys',base_resp.keys)
        if base_resp.keys != None:  # we had a response
            trials.addData('base_resp.rt', base_resp.rt)
            trials.addData('base_resp.duration', base_resp.duration)
        # the Routine "baseline" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # --- Prepare to start Routine "stim" ---
        # create an object to store info about Routine stim
        stim = data.Routine(
            name='stim',
            components=[target_stim, stim_resp],
        )
        stim.status = NOT_STARTED
        continueRoutine = True
        # update component parameters for each repeat
        # create starting attributes for stim_resp
        stim_resp.keys = []
        stim_resp.rt = []
        _stim_resp_allKeys = []
        # store start times for stim
        stim.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
        stim.tStart = globalClock.getTime(format='float')
        stim.status = STARTED
        thisExp.addData('stim.started', stim.tStart)
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
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *target_stim* updates
            
            # if target_stim is starting this frame...
            if target_stim.status == NOT_STARTED and frameN >= 0:
                # keep track of start time/frame for later
                target_stim.frameNStart = frameN  # exact frame index
                target_stim.tStart = t  # local t and not account for scr refresh
                target_stim.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(target_stim, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'target_stim.started')
                # update status
                target_stim.status = STARTED
                target_stim.setAutoDraw(True)
            
            # if target_stim is active this frame...
            if target_stim.status == STARTED:
                # update params
                pass
            
            # *stim_resp* updates
            waitOnFlip = False
            
            # if stim_resp is starting this frame...
            if stim_resp.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                stim_resp.frameNStart = frameN  # exact frame index
                stim_resp.tStart = t  # local t and not account for scr refresh
                stim_resp.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(stim_resp, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'stim_resp.started')
                # update status
                stim_resp.status = STARTED
                # keyboard checking is just starting
                waitOnFlip = True
                win.callOnFlip(stim_resp.clock.reset)  # t=0 on next screen flip
                win.callOnFlip(stim_resp.clearEvents, eventType='keyboard')  # clear events on next screen flip
            if stim_resp.status == STARTED and not waitOnFlip:
                theseKeys = stim_resp.getKeys(keyList=['up'], ignoreKeys=["escape"], waitRelease=False)
                _stim_resp_allKeys.extend(theseKeys)
                if len(_stim_resp_allKeys):
                    stim_resp.keys = _stim_resp_allKeys[-1].name  # just the last key pressed
                    stim_resp.rt = _stim_resp_allKeys[-1].rt
                    stim_resp.duration = _stim_resp_allKeys[-1].duration
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
        # check responses
        if stim_resp.keys in ['', [], None]:  # No response was made
            stim_resp.keys = None
        trials.addData('stim_resp.keys',stim_resp.keys)
        if stim_resp.keys != None:  # we had a response
            trials.addData('stim_resp.rt', stim_resp.rt)
            trials.addData('stim_resp.duration', stim_resp.duration)
        # the Routine "stim" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # --- Prepare to start Routine "delay" ---
        # create an object to store info about Routine delay
        delay = data.Routine(
            name='delay',
            components=[ISI2, delay_resp],
        )
        delay.status = NOT_STARTED
        continueRoutine = True
        # update component parameters for each repeat
        # create starting attributes for delay_resp
        delay_resp.keys = []
        delay_resp.rt = []
        _delay_resp_allKeys = []
        # store start times for delay
        delay.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
        delay.tStart = globalClock.getTime(format='float')
        delay.status = STARTED
        thisExp.addData('delay.started', delay.tStart)
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
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *delay_resp* updates
            waitOnFlip = False
            
            # if delay_resp is starting this frame...
            if delay_resp.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                delay_resp.frameNStart = frameN  # exact frame index
                delay_resp.tStart = t  # local t and not account for scr refresh
                delay_resp.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(delay_resp, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'delay_resp.started')
                # update status
                delay_resp.status = STARTED
                # keyboard checking is just starting
                waitOnFlip = True
                win.callOnFlip(delay_resp.clock.reset)  # t=0 on next screen flip
                win.callOnFlip(delay_resp.clearEvents, eventType='keyboard')  # clear events on next screen flip
            if delay_resp.status == STARTED and not waitOnFlip:
                theseKeys = delay_resp.getKeys(keyList=['up'], ignoreKeys=["escape"], waitRelease=False)
                _delay_resp_allKeys.extend(theseKeys)
                if len(_delay_resp_allKeys):
                    delay_resp.keys = _delay_resp_allKeys[-1].name  # just the last key pressed
                    delay_resp.rt = _delay_resp_allKeys[-1].rt
                    delay_resp.duration = _delay_resp_allKeys[-1].duration
                    # a response ends the routine
                    continueRoutine = False
            # *ISI2* period
            
            # if ISI2 is starting this frame...
            if ISI2.status == NOT_STARTED and frameN >= 0.0:
                # keep track of start time/frame for later
                ISI2.frameNStart = frameN  # exact frame index
                ISI2.tStart = t  # local t and not account for scr refresh
                ISI2.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(ISI2, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.addData('ISI2.started', t)
                # update status
                ISI2.status = STARTED
                ISI2.start(1*frameDur)
            elif ISI2.status == STARTED:  # one frame should pass before updating params and completing
                # Updating other components during *ISI2*
                img1.setImage(img1_file)
                img2.setImage(img2_file)
                Lcoin.setImage('input_data/coin.png')
                Rcoin.setImage('input_data/coin.png')
                Mcoin.setImage('input_data/coin.png')
                Lcross.setImage('input_data/cross.png')
                Rcross.setImage('input_data/cross.png')
                Mcross.setImage('input_data/cross.png')
                # Component updates done
                ISI2.complete()  # finish the static period
                ISI2.tStop = ISI2.tStart + 1*frameDur  # record stop time
            
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
        # check responses
        if delay_resp.keys in ['', [], None]:  # No response was made
            delay_resp.keys = None
        trials.addData('delay_resp.keys',delay_resp.keys)
        if delay_resp.keys != None:  # we had a response
            trials.addData('delay_resp.rt', delay_resp.rt)
            trials.addData('delay_resp.duration', delay_resp.duration)
        # the Routine "delay" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # --- Prepare to start Routine "task" ---
        # create an object to store info about Routine task
        task = data.Routine(
            name='task',
            components=[img1, img2, divider_line, slider_line, marker, slider_resp, task_text, submit_resp],
        )
        task.status = NOT_STARTED
        continueRoutine = True
        # update component parameters for each repeat
        # Run 'Begin Routine' code from rand_divider_slider
        leftPressed, rightPressed, marker_moved = .0, .0, .0
        positions = []
        task_txt = ''
        
        kb.clearEvents()
        divider_line.setPos([div_pos, 0])
        marker.reset()
        # create starting attributes for slider_resp
        slider_resp.keys = []
        slider_resp.rt = []
        _slider_resp_allKeys = []
        # create starting attributes for submit_resp
        submit_resp.keys = []
        submit_resp.rt = []
        _submit_resp_allKeys = []
        # store start times for task
        task.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
        task.tStart = globalClock.getTime(format='float')
        task.status = STARTED
        thisExp.addData('task.started', task.tStart)
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
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            # Run 'Each Frame' code from rand_divider_slider
            # for some reason psychopy doesnt remember this assignment
            if marker.markerPos == None:
                marker.markerPos = marker_init
            
            
            # FYI, by default press & release are false, only temporarily true
            if kb.getKeys(['left'], waitRelease=True, clear=True):
            # released
                leftPressed = 0
            if leftPressed or kb.getKeys(['left'], waitRelease=False, clear=True):
            # pressed
                leftPressed, marker_moved = 1, 1
                if (-.4 + marker_move) <= marker.markerPos:
                    marker.markerPos -= marker_move
            
            if kb.getKeys(['right'], waitRelease=True, clear=True):
            # released
                rightPressed = 0
            if rightPressed or kb.getKeys(['right'], waitRelease=False, clear=True):
            # pressed
                rightPressed, marker_moved = 1, 1
                if marker.markerPos <= (.4 - marker_move):
                    marker.markerPos += marker_move   
            
            positions.append(marker.markerPos)
            
            
            ## debug
            stimVal = 'subj_C0F1_val' if expInfo['sess_type'] in ['A','C'] else 'subj_C1F0_val'
            valence = subj_C0F1_val if stimVal == 'subj_C0F1_val' else subj_C1F0_val
            outcome, correct = 0, 0
            
            if (marker.markerPos >= div_pos and target_pos >= div_pos)\
            or (marker.markerPos <= div_pos and target_pos <= div_pos):
                
                correct = 1
                outcome = 2 if valence == 'rew' else 1
                    
            else:
                
                correct = 0
                outcome = -1 if valence == 'rew' else -2
            
            if marker_moved:
                task_txt = str(outcome)
                
                if correct and abs(target_pos - marker.markerPos) <= .05:
                    task_txt = str(outcome) + '+ bonus!' 
            
            # *img1* updates
            
            # if img1 is starting this frame...
            if img1.status == NOT_STARTED and frameN >= 0:
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
            
            # *img2* updates
            
            # if img2 is starting this frame...
            if img2.status == NOT_STARTED and frameN >= 0:
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
            
            # *divider_line* updates
            
            # if divider_line is starting this frame...
            if divider_line.status == NOT_STARTED and frameN >= 0:
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
            
            # *slider_line* updates
            
            # if slider_line is starting this frame...
            if slider_line.status == NOT_STARTED and frameN >= 0:
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
            
            # *marker* updates
            
            # if marker is starting this frame...
            if marker.status == NOT_STARTED and frameN >= 0:
                # keep track of start time/frame for later
                marker.frameNStart = frameN  # exact frame index
                marker.tStart = t  # local t and not account for scr refresh
                marker.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(marker, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'marker.started')
                # update status
                marker.status = STARTED
                marker.setAutoDraw(True)
            
            # if marker is active this frame...
            if marker.status == STARTED:
                # update params
                pass
            
            # *slider_resp* updates
            waitOnFlip = False
            
            # if slider_resp is starting this frame...
            if slider_resp.status == NOT_STARTED and frameN >= 0:
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
            if slider_resp.status == STARTED and not waitOnFlip:
                theseKeys = slider_resp.getKeys(keyList=['left','right'], ignoreKeys=["escape"], waitRelease=False)
                _slider_resp_allKeys.extend(theseKeys)
                if len(_slider_resp_allKeys):
                    slider_resp.keys = [key.name for key in _slider_resp_allKeys]  # storing all keys
                    slider_resp.rt = [key.rt for key in _slider_resp_allKeys]
                    slider_resp.duration = [key.duration for key in _slider_resp_allKeys]
            
            # *task_text* updates
            
            # if task_text is starting this frame...
            if task_text.status == NOT_STARTED and frameN >= 0:
                # keep track of start time/frame for later
                task_text.frameNStart = frameN  # exact frame index
                task_text.tStart = t  # local t and not account for scr refresh
                task_text.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(task_text, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'task_text.started')
                # update status
                task_text.status = STARTED
                task_text.setAutoDraw(True)
            
            # if task_text is active this frame...
            if task_text.status == STARTED:
                # update params
                task_text.setPos([0,.3], log=False)
                task_text.setText(task_txt, log=False)
            
            # *submit_resp* updates
            waitOnFlip = False
            
            # if submit_resp is starting this frame...
            if submit_resp.status == NOT_STARTED and frameN >= 0:
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
            if submit_resp.status == STARTED and not waitOnFlip:
                theseKeys = submit_resp.getKeys(keyList=['up'], ignoreKeys=["escape"], waitRelease=False)
                _submit_resp_allKeys.extend(theseKeys)
                if len(_submit_resp_allKeys):
                    submit_resp.keys = _submit_resp_allKeys[-1].name  # just the last key pressed
                    submit_resp.rt = _submit_resp_allKeys[-1].rt
                    submit_resp.duration = _submit_resp_allKeys[-1].duration
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
        kb.clearEvents()
        trials.addData('marker.response', marker.getRating())
        trials.addData('marker.rt', marker.getRT())
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
        # the Routine "task" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # --- Prepare to start Routine "feedback" ---
        # create an object to store info about Routine feedback
        feedback = data.Routine(
            name='feedback',
            components=[no_resp_text, Lcoin, Rcoin, Mcoin, Lcross, Rcross, Mcross, feedback_resp],
        )
        feedback.status = NOT_STARTED
        continueRoutine = True
        # update component parameters for each repeat
        # Run 'Begin Routine' code from fb_code
        # C0F1 = Curve penalty, Flat reward
        # stimVal comes from file which depends on experimenter input
        stimVal = 'subj_C0F1_val' if expInfo['sess_type'] in ['A','C'] else 'subj_C1F0_val'
        valence = subj_C0F1_val if stimVal == 'subj_C0F1_val' else subj_C1F0_val
            
        no_resp_txt = ''
        coinLR, coinM, crossLR, crossM = 0, 0, 0, 0
        correct, outcome = 0, 0
        
        if marker_moved and not submit_resp.keys == None and 'up' in submit_resp.keys:
            
            if (marker.markerPos >= div_pos and target_pos >= div_pos)\
            or (marker.markerPos <= div_pos and target_pos <= div_pos):
                
                correct = 1
                outcome = 2 if valence == 'rew' else 1
                if outcome == 2:
                    coinLR = 1
                elif outcome == 1:
                    coinM = 1
                    
            else:
                
                correct = 0
                outcome = -1 if valence == 'rew' else -2
                if outcome == -2:
                    coinLR, crossLR = 1, 1
                elif outcome == -1:
                    coinM, crossM = 1, 1   
                    
            
        elif not marker_moved:
            no_resp_txt = 'You must move the marker!'
            
        else:
            no_resp_txt = 'Respond faster!'
            
        no_resp_text.setText(no_resp_txt)
        # create starting attributes for feedback_resp
        feedback_resp.keys = []
        feedback_resp.rt = []
        _feedback_resp_allKeys = []
        # store start times for feedback
        feedback.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
        feedback.tStart = globalClock.getTime(format='float')
        feedback.status = STARTED
        thisExp.addData('feedback.started', feedback.tStart)
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
            
            # *no_resp_text* updates
            
            # if no_resp_text is starting this frame...
            if no_resp_text.status == NOT_STARTED and frameN >= 0.0:
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
            
            # *feedback_resp* updates
            waitOnFlip = False
            
            # if feedback_resp is starting this frame...
            if feedback_resp.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                feedback_resp.frameNStart = frameN  # exact frame index
                feedback_resp.tStart = t  # local t and not account for scr refresh
                feedback_resp.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(feedback_resp, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'feedback_resp.started')
                # update status
                feedback_resp.status = STARTED
                # keyboard checking is just starting
                waitOnFlip = True
                win.callOnFlip(feedback_resp.clock.reset)  # t=0 on next screen flip
                win.callOnFlip(feedback_resp.clearEvents, eventType='keyboard')  # clear events on next screen flip
            if feedback_resp.status == STARTED and not waitOnFlip:
                theseKeys = feedback_resp.getKeys(keyList=['up'], ignoreKeys=["escape"], waitRelease=False)
                _feedback_resp_allKeys.extend(theseKeys)
                if len(_feedback_resp_allKeys):
                    feedback_resp.keys = _feedback_resp_allKeys[-1].name  # just the last key pressed
                    feedback_resp.rt = _feedback_resp_allKeys[-1].rt
                    feedback_resp.duration = _feedback_resp_allKeys[-1].duration
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
        # check responses
        if feedback_resp.keys in ['', [], None]:  # No response was made
            feedback_resp.keys = None
        trials.addData('feedback_resp.keys',feedback_resp.keys)
        if feedback_resp.keys != None:  # we had a response
            trials.addData('feedback_resp.rt', feedback_resp.rt)
            trials.addData('feedback_resp.duration', feedback_resp.duration)
        # the Routine "feedback" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        thisExp.nextEntry()
        
    # completed 1.0 repeats of 'trials'
    
    if thisSession is not None:
        # if running in a Session with a Liaison client, send data up to now
        thisSession.sendExperimentData()
    
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
