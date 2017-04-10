#debug merge_line_list_altInfo
#perhaps count altInfo format has been changed

import pdb

'''
def merge_line_list_altInfo(line_list, curr_count_y_file):
    ###
    shift = 0
    if line_list[0].split()[0][0]=='[':
        shift = 1

    rel_pos = int(line_list[0].split()[shift])
    itm = line_list[0].split()[1+shift]
    res = ['%d'%rel_pos, itm]

    pdb.set_trace()

    alt_pos_list = []
    for line in line_list:
        itms = line.split()
        if len(itms)>2+shift:

            #pdb.set_trace()
            #res = res + itms[2+shift:]

            #to avoid duplicates
            candidates = itms[2+shift:]
            for c in candidates:
                gp = int(c.split(',')[0])
                #pdb.set_trace()
                if gp not in alt_pos_list:
                    #pdb.set_trace()
                    alt_pos_list.append(gp)
                    res = res + [c]
                else:
                    #pdb.set_trace()
                    continue

    #pdb.set_trace()

    res_str = '\t'.join(res)
    curr_count_y_file.write(res_str+'\n')

    return
'''

def merge_line_list_altInfo(line_list, curr_count_y_file):

    if line_list == []: return

    relPos = int(line_list[0].split()[0])
    gPos_L0 = line_list[0].split()[1]

    altInfo = {} #key: {+A,-A,+C,-C,+G,-G,+T,-T}, val: {} w/ key=alt gPos and val=lambda
    for line in line_list:
        itms = line.split()
        if len(itms)>2:
            #pdb.set_trace()
            candidates = [c.split(':') for c in itms[2].split('[') if c != ''] #e.g. candidates=[['+G', '84868432,16.0151279087,]'], ['+A', '84868432,16.0151279087,]']]
            for k, vs in candidates:#e.g. k='+G', vs='84868432,16.0151279087,]'
                if k not in altInfo: altInfo[k]={}
                gPos_L_list = vs.split(',')[:-1]#e.g. gPos_L_list=['84868432', '16.0151279087']
                for i in range(int(len(gPos_L_list)/2)):
                    gPos = int(gPos_L_list[2*i])
                    L = float(gPos_L_list[2*i+1])
                    if gPos not in altInfo[k]:
                        altInfo[k][gPos]=L
    #pdb.set_trace()

    res = []
    res.append(str(relPos))
    res.append(gPos_L0)
    #pdb.set_trace()
    st = ''
    for k, gPos_L_dic in altInfo.items():
        st += '['
        st += k + ','
        for gP, L in gPos_L_dic.items():
            st += '%d,%f,'%(gP, L)
        st += ']'
    res.append(st)
    #pdb.set_trace()
    res = '\t'.join(res)+'\n'

    curr_count_y_file.write(res)

    return res

'''
line_list = []
line_list.append('2621064 85747509,35.83  [+G:84868432,16.0151279087,][+A:84868432,16.0151279087,]')
line_list.append('2621064 85747509,35.83  [+T:84868888,1.0151279087,]')
'''

#'''
line_list = ['2324\t20642582,638285.777\t\n',
             '2324\t20642582,638285.777\t\n',
             '2324\t20642582,638285.777\t\n',
             '2324\t20642582,638285.777\t\n',
             '2324\t20642582,638285.777\t\n',
             '2324\t20642582,638285.777\t\n',
             '2324\t20642582,638285.777\t\n',
             '2324\t20642582,638285.777\t\n',
             '2324\t20642582,638285.777\t\n',
             '2324\t20642582,638285.777\t\n',
             '2324\t20642582,638285.777\t\n',
             '2324\t20642582,638285.777\t\n',
             '2324\t20642582,638285.777\t[+G:28457626,16575181.6979,23310044,37365553.6604,28902434,514919776.607,]\n',
             '2324\t20642582,638285.777\t\n',
             '2324\t20642582,638285.777\t\n',
             '2324\t20642582,638285.777\t\n',
             '2324\t20642582,638285.777\t\n',
             '2324\t20642582,638285.777\t\n',
             '2324\t20642582,638285.777\t\n',
             '2324\t20642582,638285.777\t\n',
             #'2324\t20642582,638285.777\t[+G:28457627,16575181.6980,]\n',
             ]
#'''

curr_count_y_file = open('tmp/merge_line_list_altInfo/count_y_altInfo00.txt', 'w')

print(merge_line_list_altInfo(line_list, curr_count_y_file))