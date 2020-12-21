#!/usr/bin/python
table = []
with open('ElevationModelResolution0.5.json.log') as origin_file:
    for line in origin_file:
        if line.startswith('Timer:'):
            line = line[len('Timer:'):].rstrip()
            if '--' in line:
                continue
            line = line.split()
            if line[0].startswith('dtcc'):
                table.append([line[0], 'CPU mean', 'CPU total', 'Count'])
            else:
                table.append(line)
for i in range(0,len(table)):
   if table[i][0].startswith('dtcc'):
      print(i)
   print(table[i][0],table[i][1], table[i][2], table[i][3])
