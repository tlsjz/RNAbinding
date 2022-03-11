import os
import linecache
import pickle

class dealBiochemical():
    def extractSequence(self,fastapath,pdbid):
        filelines = linecache.getlines(fastapath+'\\'+pdbid+'.fasta')
        sequenceline = filelines[1]
        sequence = sequenceline.strip('\n')
        return sequence

    def SequenceParser(self,sequence):
        encodingDic = {'H':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
                       'R':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
                       'K':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
                       'E':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
                       'D':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
                       'Y':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
                       'N':[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
                       'C':[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
                       'T':[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
                       'Q':[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
                       'S':[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
                       'W':[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
                       'M':[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
                       'F':[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
                       'P':[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                       'V':[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                       'I':[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                       'L':[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                       'A':[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                       'G':[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                       'X':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}
        hydroDic = {'H':-3.2,
                    'R':-4.5,
                    'K':-3.9,
                    'E':-3.5,
                    'D':-3.5,
                    'Y':-1.3,
                    'N':-3.5,
                    'C':2.5,
                    'T':-0.7,
                    'Q':-3.5,
                    'S':-0.8,
                    'W':-0.9,
                    'M':1.9,
                    'F':2.8,
                    'P':-1.6,
                    'V':4.2,
                    'I':4.5,
                    'L':3.8,
                    'A':1.8,
                    'G':-0.4,
                    'X':0}
        PkaDic = {'H':6,
                  'R':12.48,
                  'K':10.53,
                  'E':4.25,
                  'D':3.65,
                  'Y':10.07,
                  'N':7,
                  'C':8.18,
                  'T':7,
                  'Q':7,
                  'S':7,
                  'W':7,
                  'M':7,
                  'F':7,
                  'P':7,
                  'V':7,
                  'I':7,
                  'L':7,
                  'A':7,
                  'G':7,
                  'X':7}
        length = len(sequence)
        i = 0
        sequenceDic = {}
        while i < length:
            sequenceDic[str(i)] = []
            residuename = sequence[i]
            for each in encodingDic[residuename]:
                sequenceDic[str(i)].append(each)
            sequenceDic[str(i)].append(hydroDic[residuename])
            sequenceDic[str(i)].append(PkaDic[residuename])
            i = i+1
        return sequenceDic
        
    def process(self,windowsize,sequenceDic):
        length = len(sequenceDic.keys())
        #print(length)
        flag = 0
        gap = (windowsize-1)/2
        gap = int(gap)
        sequenceWindow = {}       
        while flag < length:
            #print(flag)
            sequenceWindow[str(flag)] = []
            if flag-gap < 0:
                for a in range(0,abs(flag-gap)):
                    for b in range(0,22):
                        sequenceWindow[str(flag)].append(0)
                c = gap - abs(flag-gap)
                for d in range(c,0,-1):
                    residuePos = str(flag-d)
                    for each in sequenceDic[residuePos]:
                        sequenceWindow[str(flag)].append(each)
                for e in range(0,gap+1):
                    residuePos = str(flag+e)
                    for each in sequenceDic[residuePos]:
                        sequenceWindow[str(flag)].append(each)
            elif length-flag <= gap:
                for f in range(gap,0,-1):
                    residuePos = str(flag-f)
                    for each in sequenceDic[residuePos]:
                        sequenceWindow[str(flag)].append(each)
                for each in sequenceDic[str(flag)]:
                    sequenceWindow[str(flag)].append(each)
                g = length-flag
                for h in range(1,g):
                    residuePos = str(flag+h)
                    for each in sequenceDic[residuePos]:
                        sequenceWindow[str(flag)].append(each)
                I = gap-(g-1)
                for j in range(0,I):
                    for k in range(0,22):
                        sequenceWindow[str(flag)].append(0)
            else:
                for l in range(gap,0,-1):
                    residuePos = str(flag-l)
                    for each in sequenceDic[residuePos]:
                        sequenceWindow[str(flag)].append(each)
                for m in range(0,gap+1):
                    residuePos = str(flag+m)
                    for each in sequenceDic[residuePos]:
                        sequenceWindow[str(flag)].append(each)
            flag = flag+1
        return sequenceWindow
    
    def run(self,windowsize,fastapath):
        filelist = os.listdir(fastapath)
        for eachfile in filelist:
            pdbid = eachfile.split('.')[0]
            sequence = self.extractSequence(fastapath,pdbid)
            sequenceDic = self.SequenceParser(sequence)
            sequenceWindow = self.process(windowsize, sequenceDic)
            picklefile = open('D:\\RNAbinding\\deep\\RB344\\sequencefeature\\'+pdbid+'.pickle','wb')
            pickle.dump(sequenceWindow,picklefile)
                    
                    
if __name__=="__main__":
    test = dealBiochemical()
    test.run(15,'D:\\RNAbinding\\RB344\\pdbfasta')
    #sequence = test.extractSequence('C:\\Users\\admin\\Desktop\\RB344\\pdbfasta','1a9n_C')
    #sequenceDic = test.SequenceParser(sequence)
    #sequenceWindow = test.process(5,sequenceDic)
    #for each in sequenceWindow.keys():
        #print(len(sequenceWindow[each]))
            
            