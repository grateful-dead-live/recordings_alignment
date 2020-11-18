from prepare import prepare_data
import json, sys
from collections import Counter
from pprint import pprint
from scipy import stats
import numpy as np
from copy import copy


DATE = '90-03-14'
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
def get_partition_bounds(points):
    parts = partition(split_segments2(points))[0]
    print('split into', len(parts))
    part_bounds = [(p[0], p[-1]) for p in parts]
    return part_bounds


def file_length(f, lengths):
    id = f.split('_')[0]
    fname = f.split('_')[1]
    return list(filter(lambda e: e[0] == fname, lengths[id]))[0][1]


def adjust_length(length, cents):
    return length / (2**(cents / 1200))


def fill_gaps(json_key, segs, jsons, lengths):
    tuning_diff = jsons[json_key]['tuning_diff']

    # interpolate gapless
    for i in range(len(segs)-1, 0, -1):
        segs.insert(i, (segs[i-1][1], segs[i][0]))

    # prepend
    first_seg_0 = segs[0][0][0]
    first_seg_1 = segs[0][0][1]
    prepended_1 = first_seg_1 - adjust_length(first_seg_0, tuning_diff)
    segs = [([0, first_seg_0], [prepended_1, first_seg_1])] + segs

    # append
    fname_0 = json_key.split('__')[0]
    last_seg_0 = segs[-1][1][0]
    last_seg_1 = segs[-1][1][1]
    l_0 = file_length(fname_0, lengths)
    l_1 = last_seg_1+adjust_length(l_0-last_seg_0, tuning_diff)
    segs += [([last_seg_0, last_seg_1], [l_0, l_1])]






    print(segs)

def main():



    subgraphs, ids_by_length, ids_by_number_of_matched_files, lengths, jsons = prepare_data(DATE)
    subgraphs = sort_subgraphs(subgraphs, lengths, ids_by_length)


    #file_length('116746_gd1990-03-14s1t02.flac', lengths)
    #sys.exit()


    #json.dump(jsons, open('jsons.json', 'w'))
    #json.dump(subgraphs, open('subgraphs.json', 'w'))
    #json.dump(lengths, open('lengths.json', 'w'))


    for n, sub in enumerate(subgraphs[1:]):
        for s in list(sub.values())[0]:
            if len(s) > 1:
                jkeys = [ track_tuple_to_json_id((s[i], s[i+1])) for i, e in  enumerate(s[:-1])]
            else:
                jkeys = [track_tuple_to_json_id((s[0], list(sub.keys())[0]))]

            dtw = jsons[jkeys[0]]['dtw']
            tuning_diff = jsons[jkeys[0]]['tuning_diff']
            #print(tuning_diff)
            file_names = jkeys[0].split('__')
            partitions = get_partition_bounds(dtw)
            print(partitions)
            #break

            #fill_gaps(jkeys[0], partitions, jsons, lengths)
            #2**(tuning_diff / 1200)
            #break

        break

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
