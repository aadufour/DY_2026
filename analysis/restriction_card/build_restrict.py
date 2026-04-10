#!/usr/bin/python
"""
  builds a restriction card by setting single operators to 1. (0.9999....)
  the ones that will not contribute to the quadratic EFT xsec will not be used for DY analysis

  usage: python3 build_restrict.py
  !! make sure to have restrict_before.txt and restrict_after.txt
  outputs all restriction cards in cards/
  the restriction cards will then have to be copied in []/mg5amcnlo/models/SMEFTsim_topU3l_MwScheme_UFO
"""



params = [
    [' 1', 'cG'],
    [' 2', 'cW'],
    [' 3', 'cH'],
    [' 4', 'cHbox'],
    [' 5', 'cHDD'],
    [' 6', 'cHG'],
    [' 7', 'cHW'],
    [' 8', 'cHB'],
    [' 9', 'cHWB'],
    [' 10', 'cuHRe'],
    [' 11', 'ctHRe'],
    [' 12', 'cdHRe'],
    [' 13', 'cbHRe'],
    [' 14', 'cuGRe'],
    [' 15', 'ctGRe'],
    [' 16', 'cuWRe'],
    [' 17', 'ctWRe'],
    [' 18', 'cuBRe'],
    [' 19', 'ctBRe'],
    [' 20', 'cdGRe'],
    [' 21', 'cbGRe'],
    [' 22', 'cdWRe'],
    [' 23', 'cbWRe'],
    [' 24', 'cdBRe'],
    [' 25', 'cbBRe'],
    [' 26', 'cHj1'],
    [' 27', 'cHQ1'],
    [' 28', 'cHj3'],
    [' 29', 'cHQ3'],
    [' 30', 'cHu'],
    [' 31', 'cHt'],
    [' 32', 'cHd'],
    [' 33', 'cHbq'],
    [' 34', 'cHudRe'],
    [' 35', 'cHtbRe'],
    [' 36', 'cjj11'],
    [' 37', 'cjj18'],
    [' 38', 'cjj31'],
    [' 39', 'cjj38'],
    [' 40', 'cQj11'],
    [' 41', 'cQj18'],
    [' 42', 'cQj31'],
    [' 43', 'cQj38'],
    [' 44', 'cQQ1'],
    [' 45', 'cQQ8'],
    [' 46', 'cuu1'],
    [' 47', 'cuu8'],
    [' 48', 'ctt'],
    [' 49', 'ctu1'],
    [' 50', 'ctu8'],
    [' 51', 'cdd1'],
    [' 52', 'cdd8'],
    [' 53', 'cbb'],
    [' 54', 'cbd1'],
    [' 55', 'cbd8'],
    [' 56', 'cud1'],
    [' 57', 'ctb1'],
    [' 58', 'ctd1'],
    [' 59', 'cbu1'],
    [' 60', 'cud8'],
    [' 61', 'ctb8'],
    [' 62', 'ctd8'],
    [' 63', 'cbu8'],
    [' 64', 'cutbd1Re'],
    [' 65', 'cutbd8Re'],
    [' 66', 'cju1'],
    [' 67', 'cQu1'],
    [' 68', 'cju8'],
    [' 69', 'cQu8'],
    [' 70', 'ctj1'],
    [' 71', 'ctj8'],
    [' 72', 'cQt1'],
    [' 73', 'cQt8'],
    [' 74', 'cjd1'],
    [' 75', 'cjd8'],
    [' 76', 'cQd1'],
    [' 77', 'cQd8'],
    [' 78', 'cbj1'],
    [' 79', 'cbj8'],
    [' 80', 'cQb1'],
    [' 81', 'cQb8'],
    [' 82', 'cjQtu1Re'],
    [' 83', 'cjQtu8Re'],
    [' 84', 'cjQbd1Re'],
    [' 85', 'cjQbd8Re'],
    [' 86', 'cjujd1Re'],
    [' 87', 'cjujd8Re'],
    [' 88', 'cjujd11Re'],
    [' 89', 'cjujd81Re'],
    [' 90', 'cQtjd1Re'],
    [' 91', 'cQtjd8Re'],
    [' 92', 'cjuQb1Re'],
    [' 93', 'cjuQb8Re'],
    [' 94', 'cQujb1Re'],
    [' 95', 'cQujb8Re'],
    [' 96', 'cjtQd1Re'],
    [' 97', 'cjtQd8Re'],
    [' 98', 'cQtQb1Re'],
    [' 99', 'cQtQb8Re'],
    [' 100', 'ceHRe'],
    [' 101', 'ceWRe'],
    [' 102', 'ceBRe'],
    [' 103', 'cHl1'],
    [' 104', 'cHl3'],
    [' 105', 'cHe'],
    [' 106', 'cll'],
    [' 107', 'cll1'],
    [' 108', 'clj1'],
    [' 109', 'clj3'],
    [' 110', 'cQl1'],
    [' 111', 'cQl3'],
    [' 112', 'cee'],
    [' 113', 'ceu'],
    [' 114', 'cte'],
    [' 115', 'ced'],
    [' 116', 'cbe'],
    [' 117', 'cje'],
    [' 118', 'cQe'],
    [' 119', 'clu'],
    [' 120', 'ctl'],
    [' 121', 'cld'],
    [' 122', 'cbl'],
    [' 123', 'cle'],
    [' 124', 'cledjRe'],
    [' 125', 'clebQRe'],
    [' 126', 'cleju1Re'],
    [' 127', 'cleQt1Re'],
    [' 128', 'cleju3Re'],
    [' 129', 'cleQt3Re'],
]

if __name__ == "__main__":

  import os
  os.makedirs ('cards', exist_ok=True)

  f_before = open ('restrict_before.txt', 'r')
  contents_before = f_before.read ()
  f_before.close ()
  f_after  = open ('restrict_after.txt', 'r')
  contents_after = f_after.read ()
  f_after.close ()

  # build restrict cards for single parameters
  # ----------------------------------------------------

  # loop over parameters to be restricted
  for param in params:
    f_restrict = open (os.path.join ('cards', 'restrict_'+ param[1] + '_massless.dat'), 'w')
    f_restrict.write (contents_before)
    
    # loop over parameters to be written
    for param2 in params:
      if param2[0] == param[0]: 
        f_restrict.write ('   ' + param2[0] + ' 9.999999e-01 # ' + param2[1] + '\n')
      else:
        f_restrict.write ('   ' + param2[0] + ' 0 # ' + param2[1] + '\n')
    # loop over parameters to be written

    f_restrict.write (contents_after)
    f_restrict.close ()
  # loop over parameters to be restricted






           