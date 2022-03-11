import os
import linecache
import pickle

class blastparser():
    def parser(self,blastoutpath):
        filelines = linecache.getlines(blastoutpath)
        length = len(filelines)
        i = 0
        seq = {}
        blastseq = {}
        querysequence = ''
        alignsequence = ''
        partnerid = ''
        score = ''
        while i < length:
            if filelines[i] == '' or filelines[i][0] != '>':
                i = i+1
            else:
                break
        #i = i+1
        #while i < length:
            #if filelines[i] == '' or filelines[i][0] != '>':
                #i = i+1
            #else:
                #break        
        #print(filelines[i])
        partnerid = filelines[i].split()[0][1:]
        #print(partnerid)
        scoreline = filelines[i+3]
        score = float(scoreline.split()[2])
        #print(score)
        contentnum = i+6
        while filelines[contentnum] != '\n':
            seqid = filelines[contentnum].split()[1]
            seq[int(seqid)] = filelines[contentnum].split()[2]
            blastid = filelines[contentnum+2].split()[1]
            blastseq[int(blastid)] = filelines[contentnum+2].split()[2]
            contentnum = contentnum+4
        #print(seq)
        #print(blastseq)
        seqIdList = []
        blastIdList = []
        seqsegment = sorted(seq.keys())
        blastsegment = sorted(blastseq.keys())
        #print(seqsegment)
        #print(blastsegment)
        segmentlen = len(seqsegment)
        segnum = 0
        while segnum < segmentlen:
            seqBegin = seqsegment[segnum]
            blastSeqBegin = blastsegment[segnum]
            SegmentSeq = seq[seqBegin]
            blastSegmentSeq = blastseq[blastSeqBegin]
            seqlength = len(SegmentSeq)
            for i in range(0,seqlength):
                if SegmentSeq[i] == blastSegmentSeq[i]:
                    if SegmentSeq[0] != '-':
                        resSeqId = seqBegin-1+i-SegmentSeq[:i].count('-')
                    else:
                        t = 0
                        while SegmentSeq[t] == '-':
                            t = t+1
                        resSeqId = seqBegin-1+i-t-SegmentSeq[t:i].count('-')
                    seqIdList.append(resSeqId)
                    if blastSegmentSeq[0] != '-':
                        blastSeqId = blastSeqBegin-1+i-blastSegmentSeq[:i].count('-')
                    else:
                        t = 0
                        while blastSegmentSeq[t] == '-':
                            t = t+1
                        blastSeqId = blastSeqBegin-1+i-t-blastsegment[t:i].count('-')
                    blastIdList.append(blastSeqId)
            segnum = segnum+1
        #print(seqIdList)
        #print(len(seqIdList))
        #print(blastIdList)
        #print(len(blastIdList))
        parserdic = {}
        parserdic['partner'] = partnerid
        parserdic['score'] = score
        parserdic['SeqId'] = seqIdList
        parserdic['BlastId'] = blastIdList
        return parserdic
    
    def siteselection(self,parserdic,sitedicpath):
        sitepickle = open(sitedicpath,'rb')
        sitedic = pickle.load(sitepickle)
        partnerId = parserdic['partner']
        pdbId = partnerId.split('_')[0]
        chainId = partnerId.split('_')[1]
        #site = sitedic[pdbId][chainId]
        site = sitedic[pdbId+'_'+chainId]
        seqIdList = parserdic['SeqId']
        blastIdList = parserdic['BlastId']
        length = len(seqIdList)
        predictedsite = []
        for i in range(0,length):
            #if str(blastIdList[i]) in site.keys():
            if blastIdList[i] in site:
                predictedsite.append(seqIdList[i])
        predict = {}
        predict['score'] = parserdic['score']
        predict['SeqId'] = predictedsite
        return predict
    
    def run(self,blastoutpath,sitedicpath):
        filelist = os.listdir(blastoutpath)
        sitepickle75 = open(sitedicpath,'rb')
        site75 = pickle.load(sitepickle75)
        for eachfile in filelist:
            pdbid = eachfile.split('.')[0]
            print(pdbid)
            #pdb = pdbid.split('_')[0]
            #chain = pdbid.split('_')[1]
            parserdic = self.parser(blastoutpath+'\\'+eachfile)
            predict = self.siteselection(parserdic, sitedicpath)
            predictpickle = open('D:\\atpbinding\\atp41\\template\\'+pdbid+'.pickle','wb')
            pickle.dump(predict,predictpickle)
            print(predict)
            #print(site75[pdbid])
if __name__=="__main__":
    test = blastparser()
    #parserdic = test.parser('D:\\atpbinding\\independentdataset\\blastout\\227out\\2XTI_B.fm0')
    #print(parserdic)
    #a = test.siteselection(parserdic,'D:\\RNAbinding\\RB344\\bindingsitedic')
    #print(a)
    #sitepickle = open('D:\\RNAbinding\\RB75\\sitedic.pickle','rb')
    #site = pickle.load(sitepickle)
    #print(site['2xzm']['B'].keys())
    #test.run('D:\\atpbinding\\atp41\\atp388blastout\\out','D:\\atpbinding\\atp388\\sitedic.pickle')
    parserdic = test.parser('D:\\atpbinding\\tool\\blastout\\templateout\\3BJU_D.fm0')
    predict = test.siteselection(parserdic, 'D:\\atpbinding\\sitedic.pickle')
    print(predict)
    #print(parserdic)
    #print(parserdic['partner'])
    #print(parserdic['BlastId'])
    #print(parserdic['SeqId'])
    

                    
            