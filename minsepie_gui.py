from gooey import Gooey, GooeyParser
import sys
import pickle
import regex as re
from Bio.SeqUtils import MeltingTemp as mt
import RNA
import os, glob
import pandas as pd
import numpy as np
import argparse
import xgboost as xgb
from pandarallel import pandarallel

from minsepie import *
from minsepie_batchinsert import *

@Gooey(program_name='MinsePIE',
      program_description = 'Predict insertion efficiencies for prime editing',
      tabbed_groups=True,
      language = 'english',
      image_dir = './img/',
      clear_before_run = True,
      terminal_font_family = 'Avenir',
      default_size=(800, 500),
    #   progress_regex=r"^Progress (\d+)$",
      menu=[{
        'name': 'About',
        'items': [{
                'type': 'AboutDialog',
                'menuTitle': 'About',
                'name': 'MinsePIE ',
                'description': 'Modelling Insertion Efficiencies for prime insertion experiments',
                'version': '1.0.2',
                'copyright': '2021',
                'website': 'https://github.com/julianeweller/MinsePIE',
                'developer': 'Juliane Weller \n Wellcome Sanger Institute \n UK',
                'license': 'Copyright (c) 2021 Juliane Weller \n Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: \n The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. \n \n THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.'
            }, {
                'type': 'Link',
                'menuTitle': 'Visit Our Lab Webite',
                'url': 'https://www.sanger.ac.uk/group/parts-group/'
            }]
          }]
      )


def ArgParse():
    """
    Asks user for all variables integrated into a GUI.
    """
    # Asks user for all variables
    # Graphical user interface
    
    parser = GooeyParser()

    subs = parser.add_subparsers(help='commands', dest='command')

    subparser_1 = subs.add_parser('SingleInsert')  
    subparser_1.add_argument('insert', metavar = 'Insert', type = str, help ='New sequences to be inserted. Do not use ambiguity characters.', gooey_options={'validator': {
                'test': "set(user_input).issubset(['A','a','T','t','G','g','C','c']) == True",
                'message': 'Only use standard DNA nucleotides: A, T, C and G.'
            }})
    subparser_1.add_argument('pbs', metavar = 'PBS', type = str, help = 'Primer binding site of pegRNA', gooey_options={'validator': {
                'test': "set(user_input).issubset(['A','a','T','t','G','g','C','c']) == True",
                'message': 'Only use standard DNA nucleotides: A, T, C and G.'
            }})
    subparser_1.add_argument('rtt', metavar = 'RTT', type = str, help = 'Reverse transcriptase template of pegRNA', gooey_options={'validator': {
                'test': "set(user_input).issubset(['A','a','T','t','G','g','C','c']) == True",
                'message': 'Only use standard DNA nucleotides: A, T, C and G.'
            }})
    subparser_1.add_argument('-m', '--mmr', metavar = 'Mismatch repair proficiency of cell line', choices=[0, 1], type = int, default = 0, help ='MMR status of cell line. 0: deficient, 1: proficient.')  #choices=[0, 1], 
    subparser_1.add_argument('-a', '--mean', metavar ='Expected mean editing rate', type = float, default = np.NAN,help ='Expected mean editing efficiency for experimental setup')
    subparser_1.add_argument('-s', '--std', metavar ='Expected standard deviation of editing rate', type = float, default = np.NAN,help ='Expected standard deviation for editing efficiency of experimental setup')

    subparser_2 = subs.add_parser('BatchMode')
    # subparser_2.add_argument('-i', '--input', dest = 'input', type = validate_table, help ='Path to csv or tsv table with insert sequences', required=True)
    subparser_2.add_argument('input', metavar ='Input file', type = str, widget = 'FileChooser', help='Path to csv or tsv table with insert sequences')
    subparser_2.add_argument('pbs', metavar = 'PBS', type = str, help = 'Primer binding site of pegRNAs')
    subparser_2.add_argument("outdir", metavar = 'Output directory', type = str, widget="DirChooser")
    subparser_2.add_argument('rtt', metavar = 'RTT', type = str, help = 'Reverse transcriptase template of pegRNAs')
    # subparser_2.add_argument('-o', '--outdir', dest = 'outdir', type = dir_path, help ='Path to output directory', required=True)
    # subparser_2.add_argument('-m', '--mmr', metavar ='Mismatch repair', type = int, default = 0, choices=[0, 1], help ='proficiency of cell line, 0: deficient, 1: proficient', widget = 'Button')
    subparser_2.add_argument('-m', '--mmr', metavar ='Mismatch repair', type = int, default = 0, choices=[0, 1], help ='proficiency of cell line, 0: deficient, 1: proficient')
    subparser_2.add_argument('-a', '--mean', metavar = 'Expected mean editing rate', type = float, default = np.NAN,help ='Expected mean editing efficiency for experimental setup')
    subparser_2.add_argument('-s', '--std', metavar = 'Expected standard deviation of editing rate', type = float, default = np.NAN,help ='Expected standard deviation for editing efficiency of experimental setup')
    
    args = parser.parse_args()


    return args

def main():
    # Initialize
    pandarallel.initialize(verbose = 1)
    args = ArgParse()

    # Load model
    packagedir = os.path.dirname(os.path.realpath(__file__))
    modeldir = packagedir + '/models/'
    try:
        assert os.path.exists(modeldir)
    except:
        raise FileNotFoundError(f'Could not find saved models. Looked in: {modeldir}')
    model_dict = load_model(modeldir)

    # Execute prediction mode
    if args.command == 'SingleInsert':
        # Create the dataframe
        request = init_df(args.insert, args.pbs, args.rtt, args.mmr, args.mean, args.std)
        request = enhance_feature_df(request)

        # Predict
        request = predict(request, model_dict)

        zscore = request['zscore'][0]
        scaledz = request['percIns_predicted'][0]
        
        if (args.mean is not np.NAN) and (args.std is not np.NAN):
            print(f'Insertion of {args.insert} \n Z-score: {zscore} \n Scaled score based on provided mean and standard deviation {scaledz}')
        else:
            print(f'Insertion of {args.insert} \n Z-score: {zscore}')
    elif args.command == 'BatchMode':
        # Bring data into shape
        if args.input.endswith('.csv'):
            request = pd.read_csv(args.input, header = None, names = ['insert'])
        elif args.input.endswith('.tsv'):
            request = pd.read_csv(args.input, header = None, sep='\t', names = ['insert'])
        else:
            print("There was an error reading in the input file. Please check the format.")
        
        request = add_batchinfo(request, args.pbs, args.rtt, args.mmr, args.mean, args.std)
        request = enhance_feature_df(request)

        # Predict
        request = predict(request, model_dict)
        
        # Display maximal inserted sequences
        n3_largest = request.nlargest(3,'zscore')

        # Save results
        outfile = os.path.basename(args.input).split('.')[0] + '_minsepie.csv'
        outpath = os.path.join(args.outdir,outfile)

        if (args.mean is not np.NAN) and (args.std is not np.NAN):
            request[['insert','zscore', 'percIns_predicted']].to_csv(outpath)
            print(n3_largest[['insert','zscore','percIns_predicted']])
        else:
            request[['insert','zscore']].to_csv(outpath)
            print(n3_largest[['insert','zscore']])
        
        print(f'Prediction results are saved as {outpath}.')
    else:
        print('Problem with mode.')

if __name__ == '__main__':
    main()