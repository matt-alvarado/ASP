"""Paths for the ASP ThermoTest vs ISORROPIA II comparison (override with environment variables)."""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
WORK  = os.environ.get('THERMO_WORK',  os.path.join(HERE, 'work'))      # run directories, ISORROPIA output, tables
OUT   = os.environ.get('THERMO_OUT',   os.path.join(HERE, 'results'))   # plots and workbook
DECKS = os.environ.get('ASP_DECKS',    os.path.join(HERE, '..', 'InorganicInputDecks'))  # default ASP deck
ISO_EXE = os.environ.get('ISO_EXE',    os.path.join(HERE, 'isobuild', 'iso'))  # ISORROPIA II built with iso_driver.f
for d in (WORK, OUT): os.makedirs(d, exist_ok=True)
