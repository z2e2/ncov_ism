import re
from datetime import date
import pandas as pd

INPUT_FOLDER = 'results'

with open('ISM_report.md') as f:
    report = f.read()

part1, part2, part3, part4 = report.split('\n<!--- dividing line --->\n')

today = date.today()
today = 'Report created on {}'.format(today.strftime("%Y/%m/%d"))
part1 = re.sub(r'Report created on [0-9|/]+', today.rstrip(), part1.rstrip())

covary_df = pd.read_csv('{}/ISM_covary_groups.txt'.format(INPUT_FOLDER))
html_covary_table = covary_df.to_html(justify='right',index=False)
covary_table = '<!--- covarying table starts --->\n{}\n<!--- covarying table ends --->'.format(html_covary_table)
part2 = re.sub(r'<!--- covarying table starts --->.*?<!--- covarying table ends --->', 
               covary_table.rstrip(), part2.rstrip(), flags=re.DOTALL)

annotation_df = pd.read_csv('{}/ISM_annotation.txt'.format(INPUT_FOLDER))
html_annotation_table = annotation_df.to_html(justify='right',index=False)
annotation_table = '<!--- annotation table starts --->\n{}\n<!--- annotation table ends --->'.format(html_annotation_table)
part3 = re.sub(r'<!--- annotation table starts --->.*?<!--- annotation table ends --->', 
               annotation_table.rstrip(), part3.rstrip(), flags=re.DOTALL)

report = '\n<!--- dividing line --->\n'.join([part1, part2, part3, part4])
with open('ISM_report.md', 'w+') as f:
    f.write(report)
