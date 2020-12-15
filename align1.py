from prepare import prepare_data
import json, sys, os
from collections import Counter
from pprint import pprint
from scipy import stats
import numpy as np
from copy import copy
from pprint import pprint
from scipy import stats
import matplotlib.pyplot as plt


#DATE = sys.argv[1]
#DATE = '90-03-14'
DATE = '90-03-24'
#DATE = '90-03-15'
MIN_R = 0.9999

flatten = lambda l: [item for sublist in l for item in sublist]

def sort_subgraphs(subs, lengths, ids_by_length):
    flatsub = lambda l: list(set([list(l.keys())[0]] + [x for sublist in list(l.values())[0] for x in sublist]))
    # file_lists sorted by list length
    file_list = []
    for i in ids_by_length:
        file_list.append([i+'_'+t[0] for t in lengths[i]])

    sorted_subs_all = []
    for i, rec in enumerate(file_list):
        sorted_subs = []
        for track in rec:
            for j, s in enumerate(subs):
                sflat = flatsub(s)
                if track in sflat:
                    sorted_subs.append(j)
                    break
        sorted_subs_all.append(list(dict.fromkeys(sorted_subs)))
    sorted_subs_all = sorted(sorted_subs_all, key=len, reverse=True)

    ordered = sorted_subs_all[0]
    for rec in sorted_subs_all[1:]:
        for i, n in enumerate(rec):
            if n not in ordered:
                prevs = [m for m in rec[:i]]
                prevs.reverse()
                nexts = [m for m in rec[i+1:]]
                pos = None
                for p in range(max([len(prevs), len(nexts)])):
                    if p < len(prevs)-1:
                        try:
                            pos = ordered.index(prevs[p]) + 1
                            break
                        except:
                            pass
                    if p < len(nexts)-1:
                        try:
                            pos = ordered.index(nexts[p])
                            break
                        except:
                            pass
                if pos != None:
                    ordered.insert(pos, n)
                else:
                    pass
                    print('cannot reorder item', rec, n)
    print(ordered)
    res = [subs[i] for i in ordered]
    return res


def track_tuple_to_json_id(n):
    return '__'.join(n)


def find_dupes(subs):
    dupes = []
    for s in subs:
        if len(list(s.values())[0]) > 0:
            dupes = [ x[0] for x in list(s.values())[0] ]
    newDict = dict(filter(lambda e: e[1] > 1, dict(Counter(dupes)).items()))
    if newDict:
        return newDict

#splits into continuous line segments
def split_segments2(points, delta=0.5):
    split = lambda l, locs: [l[i:j] for i, j in zip([0]+locs, locs+[None])]
    points = sorted(points)
    leap = lambda p, q: abs(p[0] - q[0]) > delta or abs(p[1] - q[1]) > delta
    locs = [i for i, p in enumerate(points) if i > 0 and leap(p, points[i-1])]
    return split(points, locs)

def find_best_break(segs):
    if stats.linregress(flatten(segs))[2] < MIN_R:
        rs = []
        for i in range(len(segs)):
            splits = [flatten(segs[:i]), flatten(segs[i:])]
            rs.append(max([stats.linregress(s)[2] for s in splits if len(s) > 0]))
        best_break = np.argmax(rs)
        if 0 < best_break and best_break < len(segs):
            return best_break

def partition(segs):
    partitions = [segs]
    breaks = [find_best_break(p) for p in partitions]
    while len([b for b in breaks if b is not None]) > 0:
        partitions = [[p[:breaks[i]], p[breaks[i]:]] if breaks[i] is not None else [p] for i, p in enumerate(partitions)]
        partitions = flatten(partitions)
        breaks = [find_best_break(p) for p in partitions]
        #print([len(p) for p in partitions])
    return partitions

#split into reasonably well aligned partitions
def get_partition_bounds(points, jkey):
    parts = partition(split_segments2(points))[0]
    #print(jkey, 'split into', len(parts))
    part_bounds = [[p[0], p[-1]] for p in parts]
    return part_bounds


def file_length(f, lengths):
    id = f.split('_')[0]
    fname = f.split('_')[1]  # TODO: fix: this doesn't work when underscores in filename
    return list(filter(lambda e: e[0] == fname, lengths[id]))[0][1]


def adjust_length(length, cents):
    return length / (2**(cents / 1200))

def plotFigure(segs, json_key, lengths, fname, dtw, jsons):
    
    p = plt.figure()
    for s in segs:
        #print(s[0], s[1])
        try:
            colour = s[2]
            plt.plot([s[0][0], s[1][0]], [s[0][1], s[1][1]], color=colour, alpha=0.5)
        except:
            end = dtw.index(s[1]) + 1
            start = dtw.index(s[0])
            sdtw = np.array(dtw)[start:end]
            #plt.plot(sdtw[:,0], sdtw[:,1], color='b', alpha=0.5)  # plot all values from dtw between start and end
            plt.plot([s[0][0], s[1][0]], [s[0][1], s[1][1]], color='b', alpha=0.5) # plot line from start to end

            #x = sdtw[:,0]
            #y = sdtw[:,1]
            #coef = np.polyfit(x,y,1)
            #poly1d_fn = np.poly1d(coef)
            #plt.plot(x, poly1d_fn(x), color='b', alpha=0.5)    # plot linear regression line of dtw segment

            #ratio =  1 / 2**(jsons[json_key]['tuning_diff'] / 1200)
            #len1 = s[1][0] - s[0][0]
            #print(len1)
            #print(ratio)


            #plt.plot([s[0][0], s[1][0]], [s[0][1], s[0][1] + len1 * ratio], color='b', alpha=0.5)
            #sys.exit()


        '''
        # use all points of dtw?
        if colour == 'r':
            pass
            #plt.plot([s[0][0], s[1][0]], [s[0][1], s[1][1]], color=colour, alpha=0.5)
        else:
            points = dtw[dtw.index([s[0][0], s[1][0]]):dtw.index([s[0][1], s[1][1]])]
            print(points)
        '''

        #break
    plt.tight_layout()

    p.savefig(fname+'.pdf', bbox_inches='tight')
    plt.close(p)



def fill_gaps(json_key, segs, jsons, lengths):
    tuning_diff = jsons[json_key]['tuning_diff'] 
    #print(segs)
    new_segs = copy(segs)
    for n, s in enumerate(segs):
        if n == 0:
            start = 0
        else:
            start = segs[n-1][1][0]
        end = s[0][0]
        pre_seg = [[[start, s[0][1]-adjust_length(end-start, tuning_diff)], [s[0][0], s[0][1]], 'r']]
        new_segs = pre_seg + new_segs
        #print(n, pre_seg)

        if n == len(segs)-1:
            fname_0 = json_key.split('__')[0] 
            l_0 = file_length(fname_0, lengths)
            #print(l_0)
            app_seg = [[[s[1][0], s[1][1]], [l_0, s[1][1]+adjust_length(l_0-s[1][0], tuning_diff)], 'r' ]]
            new_segs += app_seg
            #print(n, app_seg)
    return sorted(new_segs, key=lambda s: s[0][0])

    
'''
def linReg(partitions):
    print(partitions)
    p = [[i for i in j] for j in partitions[0]]

    init_slope, intercept, init_r_value, p_value, std_err = stats.linregress(np.swapaxes(p[0],0,1))
    print('first', init_slope, init_r_value**2)

    for n, i in enumerate(p[:-1]):
        q = np.swapaxes(np.vstack(p[n:n+2]),0,1)
        slope, intercept, r_value, p_value, std_err = stats.linregress(q)
        #print(q)
        print(n, slope, r_value**2)

    #test: 2nd segment with 10s jump:
    q[1][2:] += 10
    slope, intercept, r_value, p_value, std_err = stats.linregress(q)
    #print(q)
    print('test', slope, r_value**2)
'''

'''
def process_chain_BAK(c, all_partitions, partition_jkeys):
    
    #pprint(partition_jkeys)
    #return
    translation = []
    # try for length 2 only:
    for i, t in enumerate(c[1][:-1]):
        jk = track_tuple_to_json_id((t, c[1][i+1])) 
        translation.append(all_partitions[partition_jkeys.index(jk)])
    json.dump(translation, open('translation.json', 'w'))



    first_seg = translation[0][2][:2]
    print()
    print('original:  ', first_seg)
    inter = translation[1]

    match_start_seg = list(filter(lambda x: x[0][0] <= first_seg[0][1] <= x[1][0], inter))[0]#[:2]
    match_end_seg = list(filter(lambda x: x[0][0] <= first_seg[1][1] <= x[1][0], inter))[0]#[:2]
    
    start_chain = [match_start_seg]
    end_chain = [match_end_seg]

    def map_seg(p, s):
        prop = (p - s[0][0]) / (s[1][0] - s[0][0])
        return prop * (s[1][1] - s[0][1]) + s[0][1]
    
    def map_chain(seg, chain):
        for s in chain:
            seg = [seg[0], map_seg(seg[1], s)]
        return seg

    start_to_ref = map_chain(first_seg[0], start_chain)
    end_to_ref = map_chain(first_seg[1], end_chain)

    print('start seg: ', match_start_seg)
    print('end seg:   ', match_end_seg)
    print('new seg:   ', start_to_ref, end_to_ref)
'''



def process_chain(c, all_partitions, partition_jkeys, jsons):
    def map_seg(p, s):
            prop = (p - s[0][0]) / (s[1][0] - s[0][0])
            return prop * (s[1][1] - s[0][1]) + s[0][1]
    
    print()
    print()
    print(c)
    jk1 = track_tuple_to_json_id((c[0], c[1])) 
    jk2 = track_tuple_to_json_id((c[1], c[-1])) 
    translation = [all_partitions[partition_jkeys.index(jk1)], all_partitions[partition_jkeys.index(jk2)]]

    #tuning_diff = jsons[jk1]['tuning_diff'] + jsons[jk2]['tuning_diff']


    new_segments = []
    for s in translation[0]:
        pend_flag = False
        seg = s[:2]
        print('original:  ', seg)


        # if start[0][1] < 0 or start[1][1] > length the first/last segment is used. tuning_diff instead?
        search_segment = list(filter(lambda x: x[0][0] <= seg[0][1] <= x[1][0], translation[1]))
        if search_segment:
            match_start_seg = search_segment[0][:2]
            print('start seg: ', match_start_seg)
        else:
            print('prepend to ', translation[1][0][:2])
            match_start_seg = translation[1][0][:2]
            #match_start_seg = [[0, adjust_length(1, jsons[jk2]['tuning_diff'])+translation[1][0][0][1]], translation[1][0][1]]  # ??
            print('prepend to ', match_start_seg)
            pend_flag = True
     
        search_segment = list(filter(lambda x: x[0][0] <= seg[1][1] <= x[1][0], translation[1]))
        if search_segment:
            match_end_seg = search_segment[-1][:2]
            print('end seg:   ', match_end_seg)
        else:
            print('append to  ', translation[1][-1][:2])
            end_to_ref = translation[1][-1][:2]
            pend_flag = True
           
        
        start_to_ref = [seg[0][0], map_seg(seg[0][1], match_start_seg)]
        end_to_ref = [seg[1][0], map_seg(seg[1][1], match_end_seg)]

        new_segment = [start_to_ref, end_to_ref, 'c']
        new_segments.append(new_segment)
        
        print('new seg:   ', new_segment)
        print()
        #sys.exit()

    #  TODO: fill gap for new partition
    all_partitions.append(new_segments)
    new_jkey = track_tuple_to_json_id((c[0], c[-1]))
    partition_jkeys.append(new_jkey)

    return all_partitions, partition_jkeys
    



def main():

    subgraphs, ids_by_length, ids_by_number_of_matched_files, lengths, jsons = prepare_data(DATE)
    subgraphs = sort_subgraphs(subgraphs, lengths, ids_by_length)


    #file_length('116746_gd1990-03-14s1t02.flac', lengths)
 


    #json.dump(jsons, open('jsons.json', 'w'))
    #json.dump(lengths, open('lengths.json', 'w'))
    json.dump(subgraphs, open('subgraphs.json', 'w'))
    #sys.exit()
    all_partitions = []
    partition_jkeys = []
    for n, sub in enumerate(subgraphs[16:]):
        chains = [] # json keys of chained alignments
        for s in list(sub.values())[0]:
            #if len(s) > 1:
            if len(s) > 1:
                jkeys = [ track_tuple_to_json_id((s[i], s[i+1])) for i, e in enumerate(s[:-1])]
                #chains.append((len(s), s + list(sub.keys())))
                chains.append(s + list(sub.keys()))
            else:
                jkeys = [track_tuple_to_json_id((s[0], list(sub.keys())[0]))]
            #print(jkeys[0])
            dtw = jsons[jkeys[0]]['dtw']
            dtw = [[x[1], x[0]] for x in dtw] #swap columns to match order of file names/lengths
            tuning_diff = jsons[jkeys[0]]['tuning_diff']
            #print(tuning_diff)
            #file_names = jkeys[0].split('__')
            partitions = get_partition_bounds(dtw, jkeys[0])
            
            #all_partitions.append(fill_gaps(jkeys[0], partitions, jsons, lengths))
            #2**(tuning_diff / 1200)
            partitions = fill_gaps(jkeys[0], partitions, jsons, lengths)
            #print(partitions)
            all_partitions.append(partitions)
            partition_jkeys.append(jkeys[0])
            #with open('plots/'+jkeys[0]+'.txt', 'w') as sfile:
            #    pprint(partitions, sfile)

            target_folder = os.path.join('plots', DATE)
            if not os.path.exists(target_folder):
                os.mkdir(target_folder)

            fname = f'{target_folder}/{jkeys[0]}'
            #print(fname)
            #json.dump(sorted(partitions, key=lambda x: x[0][0]), open(fname+'.json', 'w'))
            #sys.exit()
            #plotFigure(partitions, jkeys[0], lengths, fname, dtw, jsons)

            #break
    
        for c in chains:
            all_partitions, partition_jkeys = process_chain(c, all_partitions, partition_jkeys, jsons)
            #break
        #json.dump(all_partitions, open('all_partition.json', 'w'))
        break
    json.dump(all_partitions, open('all_partition.json', 'w'))

    
    '''
    #find overlaps in reference
    for j, p in enumerate(all_partitions):
        for i in range(len(p)):
            if i > 0 and p[i-1][1][1] > p[i][0][1]+1: #starts earlier on the y axis
                #there's an overlap
                print(partition_jkeys[j])
                print(p[i-1], p[i])
                print()
    '''          

        #all_partitions = linReg(all_partitions)
            
            #with open(jkeys[0] + '.txt', 'w') as sfile:
            #.l    pprint(all_partitions[0], sfile)
            #plotFigure(all_partitions[0], jkeys[0], lengths)
            
            
    
    
    #d = find_dupes(subgraphs)
    #pprint('dupes:', d)


main()


'''
{
  "116746_gd1990-03-14s1t02.flac": [
    ["116030_GD90-03-14d1t02.flac"],
    ["89689_gd1990-03-14d1t01.flac"],
    ["83778_gd1990-03-14d1t02.flac", "116030_GD90-03-14d1t02.flac"],
    ["125852_gd1990-03-14.Nak300.t01.flac", "89689_gd1990-03-14d1t01.flac"]
  ]
}
'''
