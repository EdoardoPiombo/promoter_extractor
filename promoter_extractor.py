import argparse as ap
import pandas as pd
from Bio import SeqIO

#arguments
parser = ap.ArgumentParser()
parser.add_argument("-g", "--gff", help="input gff file")
parser.add_argument("-s", "--genome", help="input genome fasta file")
parser.add_argument("-l", "--length", help="maximum promoter length")
parser.add_argument("-p", "--prevent_overlapping",action='store_true', help="prevents promoters of different genes to overlap one another. When two promoters would overlap each other the distance between them is equally distributed")
parser.add_argument("-o", "--output", help="output gff file")

args = parser.parse_args()


gff = args.gff
genome = args.genome
output = args.output
length = int(args.length)

data = pd.read_csv(gff, sep = '\t', names = ['scaffold', 'predictor', 'feature', 'start', 'end', 'thing', 'sense', 'RF', 'description'], comment = '#')
data = data.where(data['feature'] == 'gene').dropna()
scaflist = data['scaffold'].drop_duplicates().tolist()
data['to_remove'] = 'No'

# 2 dataframes are created, one for promoters and one for genes, and then they are both used to create the final one
finaldata = pd.DataFrame(columns = data.columns)

finalid = -1
for record in SeqIO.parse(open(genome), "fasta"):
    scafdata = data.where(data['scaffold'] == record.id).dropna()
    scafdata = scafdata.reset_index(drop = True)
    scafpromdata = scafdata.copy()
    scafpromdata['feature'] = 'promoter'
    for id in scafdata.index.tolist():
        promline = scafpromdata.loc[id]
        line = scafdata.loc[id]
        if line['sense'] == '+':
            if line['start'] > length :
                promline['end'] = promline['start'] -1
                promline['start'] -= length
                line['start'] = promline['start']
            elif line['start'] <= length :
                promline['end'] = promline['start'] -1
                promline['start'] = 1
                # genes that start at the first base have no promoter as there are no bases before the gene
                if line['start'] == 1:
                    line['to_remove'] = 'yes'
                    promline['to_remove'] = 'yes'
                line['start'] = promline['start']
            if int(id) > 0:
                prevline = scafdata.loc[id-1]
                promprevline = scafpromdata.loc[id-1]
                if prevline['sense'] == '+' :
                    if promline['start'] <= prevline['end'] :
                        promline['start'] = prevline['end'] + 1
                        line['start'] = promline['start']
                        if promline['start'] > promline['end']:
                            promline['start'] = promline['end'] - 1
                            line['start'] = promline['start']
                elif prevline['sense'] == '-':
                    try:
                        if args.prevent_overlapping == True:
                    # This part is to make sure to not have overlapping promoters.
                    # When two promoters meet each other the distance between them is equally distributed.
                    # The finaldata of the previous line is changed accordingly
                            if promline['start'] <= prevline['end'] :
                                distance = promline['end'] - promprevline['start']
                                if distance > 0:
                                    promprevline['end'] = promprevline['start'] + (round(distance/2))
                                    prevline['end'] = promprevline['end']
                                    finaldata.loc[finalid-1] = promprevline
                                    finaldata.loc[finalid] = prevline
                                    promline['start'] = promline['end'] - (round(distance/2)) + 1
                                    line['start'] = promline['start']
                                    if promline['start'] > promline['end']:
                                        promline['start'] = promline['end'] - 1
                                        line['start'] = promline['start']
                                else:
                                    prevline['end'] = promprevline['start'] + 1
                                    promprevline['end'] = prevline['end']
                                    promline['start'] = promline['end'] - 1
                                    line['start'] = promline['start']
                    except AttributeError:
                        pass
        elif line['sense'] == '-' :
            if (line['end'] + length) < len(record.seq):
                promline['start'] = promline['end'] +1
                promline['end'] += length
                line['end'] = promline['end']
            else:
                promline['start'] = promline['end'] +1
                promline['end'] = len(record.seq)
                # genes on the minus filament which start at the end of the scaffold have no promoter as there are no bases after the gene
                if line['end'] == len(record.seq) :
                    line['to_remove'] = 'yes'
                    promline['to_remove'] = 'yes'
                line['end'] = promline['end']                
            try:
                nextline = scafdata.loc[id+1]
                if nextline['sense'] == '+':
                    if line['end'] > nextline['start']:
                        line['end'] = nextline['start'] -1
                        promline['end'] = line['end']
                        if promline['start'] > promline['end']:
                            promline['end'] = promline['start'] + 1
                            line['end'] = promline['end']
                elif nextline['sense'] == '-':
                    if line['end'] >= nextline['start']:
                        line['end'] = nextline['start'] -1
                        promline['end'] = nextline['start'] -1
                        if promline['start'] > promline['end']:
                            promline['end'] = promline['start'] + 1
                            line['end'] = promline['end']
            except KeyError:
                pass
        if line['start'] <= 0:
            line['start'] = 1
        if promline['start'] <= 0:
            promline['start'] = 1
        if promline['start'] >= promline['end']:
            promline['end'] = promline['start']
        scafdata.loc[id] = line
        scafpromdata.loc[id] = promline
        finalid +=2
        finaldata.loc[finalid-1] = promline
        finaldata.loc[finalid] = line

droppable = finaldata.where(finaldata['to_remove'] == 'yes').dropna().index.to_list()
finaldata.drop(['to_remove'], axis = 1, inplace = True)
finaldata.drop(droppable, axis = 0, inplace = True)
finaldata = finaldata.reset_index(drop = True)
        
finaldata.to_csv(path_or_buf=output, sep = '\t', index =False, header=False)
