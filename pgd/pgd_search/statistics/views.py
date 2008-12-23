from django.shortcuts import render_to_response
from django.template import RequestContext
from django.conf import settings
from statlib import stats

from constants import AA_CHOICES
from pgd_search.models import Search, Segment, iIndex

"""
display statistics about the search
"""
def searchStatistics(request):

    stat_attributes_base = ['L1','L2','L3','L4','L5','a1','a2','a3','a4','a5','a6','a7']
    TOTAL_INDEX = {'na':0,'e':1,'E':2,'S':3,'h':4,'H':5,'t':6,'T':7,'g':8,'G':9,'B':10,'i':11,'I':12}
    STAT_INDEX = {}

    ss_field = 'r%i_ss' % iIndex
    aa_field = 'r%i_aa' % iIndex

    stat_attributes = ['r%i_%s'%(iIndex, f) for f in stat_attributes_base]
    fieldNames      = ['r%i_%s'%(iIndex, f) for f in stat_attributes_base]
    fieldNames.append(ss_field)

    for i in range(len(stat_attributes)):
        STAT_INDEX[stat_attributes[i]] = i

    # get search from session
    search = request.session['search']
    searchQuery = search.querySet()

    peptides = {}
    #iterate through the aa_choices
    for code,long_code in AA_CHOICES:
        #create data structure
        peptide = {
                'longCode':long_code,
                #total
                'total':0,
                # attributes with just sums
                'counts':[['na',0],['e',0],['E',0],['S',0],['h',0],['H',0],['t',0],['T',0],['g',0],['G',0],['B',0],['i',0],['I',0]],
                # attributes with stats
                'stats':[['L1',[]],['L2',[]],['L3',[]],['L4',[]],['L5',[]],['a1',[]],['a2',[]],['a3',[]],['a4',[]],['a5',[]],['a6',[]],['a7',[]]]
            }

        #query segments matching this AA with just the fields we want to perform calcuations on
        residueData = searchQuery.filter(**{aa_field:code}).values(*fieldNames)
        peptide['total'] = searchQuery.filter(**{aa_field:code}).count()

        #iterate through all the segment data
        for data in residueData:

            #calculate values
            if data[ss_field] == ' ':
                peptide['counts'][TOTAL_INDEX['na']][1] += 1
            else:
                peptide['counts'][TOTAL_INDEX[data[ss_field]]][1] += 1

            #store all values for attributes into arrays
            for key in stat_attributes:
                peptide['stats'][STAT_INDEX[key]][1].append(data[key])

        #calculate statistics
        for attribute in stat_attributes:
            list = peptide['stats'][STAT_INDEX[attribute]][1]
            list_len = len(list)
            if list_len > 1:
                mean = stats.mean(list)
                #now that we have mean calculate standard deviation
                stdev = stats.stdev(list)
                range_min = '%+.3f' % (min(list) - mean)
                range_max = '%+.3f' % (max(list) - mean)

            # if theres only 1 item then the stats are simpler to calculate
            elif list_len == 1:
                mean = list[0]
                std_dev = 0
                range_min = 0
                range_max = 0

            else:
                mean = 0
                stdev = 0
                range_min = 0
                range_max = 0

            peptide['stats'][STAT_INDEX[attribute]][1] = {'mean':mean,'std':stdev,'min':range_min, 'max':range_max}

        peptides[code] = peptide

    return render_to_response('stats.html', {
        'attributes': stat_attributes_base,
        'peptides':peptides
    }, context_instance=RequestContext(request))