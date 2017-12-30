from requests import ConnectionError
import logging
import subprocess
try:
  import json
  import gspread
  from oauth2client.client import SignedJwtAssertionCredentials
except ImportError:
  logging.warning("Results cannot be uploaded on the website - Module missing")
import os
import time

from Paths import *
import Parameters


def upload(cands,idL):
  upload_sheet(cands,idL)
  upload_plots(idL)
  return


def upload_plots(idL):
  folder = '{}/sp/candidates'.format(TEMP_FOLDER.format(idL))
  FNULL = open(os.devnull, 'w')
  error = subprocess.call(['scp','-prq',folder,'ag004:/var/www/lofarpwg/lotaas-sp/observations/{}/'.format(idL)], stdout=FNULL, stderr=FNULL)
  if error: raise ConnectionError("ATTENTION! Website currently down. Try to upload the observation later with Upload.py")
  return


def upload_sheet(cands,idL):
  try: json_key = json.load(open(SITE_CERT))
  except IOError:
    logging.warning("Spreadsheet cannot be uploaded - Google certificate missing")
    return
  scope = ['https://spreadsheets.google.com/feeds']
  credentials = SignedJwtAssertionCredentials(json_key['client_email'], json_key['private_key'], scope)
  gc = gspread.authorize(credentials)
  sh = gc.open("LSP_candidates")

  wks = sh.worksheet("Candidates")
  col_size = wks.col_count
  col_size_letter = chr(col_size-1 + ord('A'))
  row_size = wks.row_count

  #Create a new backup sheet
  '''
  backup = sh.worksheet("Backup")
  sh.del_worksheet(backup)
  backup = sh.add_worksheet(title="Backup", rows=row_size, cols=col_size)
  old_cells = backup.range('A1:{col}{row}'.format(row=row_size,col=col_size_letter))
  new_cells = wks.range('A1:{col}{row}'.format(row=row_size,col=col_size_letter))
  for i in range(len(old_cells)):         
    old_cells[i].value = new_cells[i].input_value
  backup.update_cells(old_cells)
  '''

  #Remove old rows for the same observation
  obs_list = wks.col_values(4)
  row_list = [i+1 for i,val in enumerate(obs_list) if val == idL]
  for row_num in row_list:
    cells = wks.range('A{row}:{col}{row}'.format(row=row_num,col=col_size_letter))
    for cell in cells : cell.value = ''
    wks.update_cells(cells)

  #Append candidates
  date = time.strftime("%m/%d/%Y")
  fileName, fileExtension = os.path.splitext(Parameters.__file__)
  git_folder = '{}/.git'.format(os.path.dirname(os.path.dirname(fileName)))
  vers = subprocess.check_output(['git','--git-dir',git_folder,'describe','--tags','--abbrev=0','--always']).strip()
  for idx,cand in cands.iterrows():
    if cand.N_pulses == 1: kind = 'SP'
    else: kind = 'RC'
    link = '=HYPERLINK(CONCATENATE("http://www.astron.nl/lofarpwg/lotaas-sp/observations/{}/";OFFSET($A$1;ROW()-1;0);".png");"Plot")'.format(idL)
    row = [cand.id, date, vers, idL, cand.SAP, cand.BEAM, kind, cand.N_pulses, cand.DM, cand.Rank, '', '', 'ToProcess', '', link]
    wks.append_row(row)

  #Sort spreadsheet
  wks = sh.worksheet("Candidates")
  col_size = wks.col_count
  col_size_letter = chr(col_size-1 + ord('A'))
  row_size = wks.row_count
  old_cells = wks.range('A1:{col}{row}'.format(row=row_size,col=col_size_letter))
  new_cells = list(old_cells)
  old_col = wks.col_values(1)
  new_idx = [old_col.index(i) for i in sorted(old_col, reverse=True)]
  for i, i_new in enumerate(new_idx):
    for col in range(col_size):
      new_cells[i*col_size+col].value = old_cells[i_new*col_size+col].input_value
  wks.update_cells(new_cells)

  #Resize spreadsheet
  try:
    row = wks.col_values(1)
    row = [val for val in row if val != '']
    row = len(row)
    wks.resize(rows=row, cols=col_size)
  except ValueError: pass

  return
