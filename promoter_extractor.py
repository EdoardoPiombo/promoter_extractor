import argparse as ap
import pandas as pd
from Bio import SeqIO

#arguments
parser = ap.ArgumentParser()
parser.add_argument("-g", "--gff", help="input gff file")
parser.add_argument("-d", "--genome", help="input genome fasta file")
parser.add_argument("-d", "--length", help="maximum promoter length")
parser.add_argument("-o", "--output", help="output gff file")

args = parser.parse_args()


gff = args.gff
genome = args.genome
output = args.output
length = args.length

data = pd.read_csv(gff, sep = '\t', names = ['scaffold', 'predictor', 'feature', 'start', 'end', 'thing', 'sense', 'RF', 'description'])
data = data.where(data['feature'] == 'gene').dropna()
scaflist = data['scaffold'].drop_duplicates().tolist()

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
        elif line['sense'] == '-' :
            if (line['end'] + length) < len(record.seq):
                promline['start'] = promline['end'] +1
                promline['end'] += length
                line['end'] = promline['end']
            else:
                promline['start'] = promline['end'] +1
                promline['end'] = len(record.seq)
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
            promline['end'] = promline['start'] + 1
        scafdata.loc[id] = line
        scafpromdata.loc[id] = promline
        finalid +=2
        finaldata.loc[finalid-1] = promline
        finaldata.loc[finalid] = line
        
finaldata = finaldata.reset_index(drop = True)
        
finaldata.to_csv(path_or_buf=output, sep = '\t', index =False, header=False)