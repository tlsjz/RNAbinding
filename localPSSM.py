import math
import linecache
import os
import pickle


class dealPSSM():
    def PSSMparser(self,pssmpath,pdbid):
        filelines = linecache.getlines(pssmpath+'\\'+pdbid+'.pssm')
        pssmDic = {}
        for line in filelines:
            content = line.split()
            if len(content) == 44:
                residuePosition = int(content[0])-1
                pssmDic[str(residuePosition)] = []
                for i in range(2,22):
                    #pssmDic[str(residuePosition)].append(int(content[i]))
                    pssmDic[str(residuePosition)].append(self.normalize(int(content[i])))
        return pssmDic
    
    def normalize(self,value):
        a = 1+math.exp(value)
        b = 1/a
        return b
    
    def windowPSSM(self,windowsize,PSSMdic):
        length = len(PSSMdic.keys())
        gap = (windowsize-1)/2
        gap = int(gap)
        flag = 0
        PSSMwindow = {}        
        while flag < length:
            PSSMwindow[str(flag)] = []
            if flag-gap < 0:
                for a in range(0,abs(flag-gap)):
                    for b in range(0,20):
                        PSSMwindow[str(flag)].append(0)
                c = gap - abs(flag-gap)
                for d in range(c,0,-1):
                    residuePos = str(flag-d)
                    for each in PSSMdic[residuePos]:
                        PSSMwindow[str(flag)].append(each)
                for e in range(0,gap+1):
                    residuePos = str(flag+e)
                    for each in PSSMdic[residuePos]:
                        PSSMwindow[str(flag)].append(each)
            elif length-flag <= gap:
                for f in range(gap,0,-1):
                    residuePos = str(flag-f)
                    for each in PSSMdic[residuePos]:
                        PSSMwindow[str(flag)].append(each)
                for each in PSSMdic[str(flag)]:
                    PSSMwindow[str(flag)].append(each)
                g = length-flag
                for h in range(1,g):
                    residuePos = str(flag+h)
                    for each in PSSMdic[residuePos]:
                        PSSMwindow[str(flag)].append(each)
                I = gap-(g-1)
                for j in range(0,I):
                    for k in range(0,20):
                        PSSMwindow[str(flag)].append(0)
            else:
                for l in range(gap,0,-1):
                    residuePos = str(flag-l)
                    for each in PSSMdic[residuePos]:
                        PSSMwindow[str(flag)].append(each)
                for m in range(0,gap+1):
                    residuePos = str(flag+m)
                    for each in PSSMdic[residuePos]:
                        PSSMwindow[str(flag)].append(each)         
            flag = flag+1
        return PSSMwindow
    
    def pairRelation(self,windowsize,PSSMdic):
        length = len(PSSMdic.keys())
        gap = (windowsize-1)/2
        gap = int(gap)
        flag = 0  
        pairList = {}
        while flag < length:
            pairList[str(flag)] = []
            for a in range(0,400):
                pairList[str(flag)].append(0)            
            if flag-gap < 0:
                b = gap - abs(flag-gap)
                for c in range(b,0,-1):
                    residuePos = str(flag-c)
                    tempresult = self.listmultiply(PSSMdic[str(flag)],PSSMdic[str(residuePos)])
                    pairList[str(flag)] = self.listplus(pairList[str(flag)],tempresult)
                for d in range(1,gap+1):
                    residuePos = str(flag+d)
                    tempresult = self.listmultiply(PSSMdic[str(flag)],PSSMdic[str(residuePos)])
                    pairList[str(flag)] = self.listplus(pairList[str(flag)],tempresult)
            elif length-flag <= gap:
                for e in range(gap,0,-1):
                    residuePos = str(flag-e)
                    tempresult = self.listmultiply(PSSMdic[str(flag)],PSSMdic[str(residuePos)])
                    pairList[str(flag)] = self.listplus(pairList[str(flag)],tempresult)
                f = length-flag
                for g in range(1,f):
                    residuePos = str(flag+g)
                    tempresult = self.listmultiply(PSSMdic[str(flag)],PSSMdic[str(residuePos)])
                    pairList[str(flag)] = self.listplus(pairList[str(flag)],tempresult)
            else:
                for h in range(gap,0,-1):
                    residuePos = str(flag-h)
                    tempresult = self.listmultiply(PSSMdic[str(flag)],PSSMdic[str(residuePos)])
                    pairList[str(flag)] = self.listplus(pairList[str(flag)],tempresult)
                for i in range(1,gap+1):
                    residuePos = str(flag+i)
                    tempresult = self.listmultiply(PSSMdic[str(flag)],PSSMdic[str(residuePos)])
                    pairList[str(flag)] = self.listplus(pairList[str(flag)],tempresult)
            flag = flag+1
        return pairList
                
    def listmultiply(self,a,b):
        c = []
        for eacha in a:
            for eachb in b:
                c.append(eacha*eachb)
        return c   
    
    def listplus(self,a,b):
        c = []
        for i in range(len(a)):
            c.append(a[i]+b[i])
        return c    
    
    def multiRelation(self,windowsize,PSSMdic):
        length = len(PSSMdic.keys())
        gap = (windowsize-1)/2
        gap = int(gap)
        flag = 0  
        multiDic = {}  
        while flag < length:
            leftList = []
            rightList = []
            for a in range(0,20):
                leftList.append(0)
                rightList.append(0)
            if flag-gap < 0:
                b = gap - abs(flag-gap)
                for c in range(b,0,-1):
                    residuePos = str(flag-c)
                    leftList = self.listplus(leftList,PSSMdic[str(residuePos)])
                leftList = self.listplus(leftList,PSSMdic[str(flag)])
                for d in range(0,gap+1):
                    residuePos = str(flag+d)
                    rightList = self.listplus(rightList,PSSMdic[str(residuePos)])
            elif length-flag <= gap:
                for e in range(gap,0,-1):
                    residuePos = str(flag-e)
                    leftList = self.listplus(leftList,PSSMdic[str(residuePos)])
                leftList = self.listplus(leftList,PSSMdic[str(flag)])
                f = length-flag
                for g in range(0,f):
                    residuePos = str(flag+g)
                    rightList = self.listplus(rightList,PSSMdic[str(residuePos)])
            else:
                for h in range(gap,0,-1):
                    residuePos = str(flag-h)
                    leftList = self.listplus(leftList,PSSMdic[str(residuePos)])
                leftList = self.listplus(leftList,PSSMdic[str(flag)])
                for i in range(1,gap+1):
                    residuePos = str(flag+i)
                    rightList = self.listplus(rightList,PSSMdic[str(flag)])
            multiDic[str(flag)] = leftList+rightList
            flag = flag+1
        return multiDic
    
    def run(self,fastapath,pssmpath,windowsize):
        filelist = os.listdir(fastapath)
        for eachfile in filelist:
            pssmfeature = {}
            pdbid = eachfile.split('.')[0]
            PSSMdic = self.PSSMparser(pssmpath,pdbid)
            pssmwindow = self.windowPSSM(windowsize,PSSMdic)
            pairdic = self.pairRelation(windowsize, PSSMdic)
            multidic = self.multiRelation(windowsize, PSSMdic)
            length = len(PSSMdic.keys())
            for i in range(0,length):
                pssmfeature[str(i)] = []
                for each in pssmwindow[str(i)]:
                    pssmfeature[str(i)].append(each)
                for each in pairdic[str(i)]:
                    pssmfeature[str(i)].append(each)
                for each in multidic[str(i)]:
                    pssmfeature[str(i)].append(each)
            picklefile = open('D:\\RNAbinding\\deep\\RB344\\pssmfeature\\'+pdbid+'.pickle','wb')
            pickle.dump(pssmfeature,picklefile)   
                         
if __name__=="__main__":
    test = dealPSSM()
    test.run('D:\\RNAbinding\\RB344\\pdbfasta','D:\\RNAbinding\\RB344\\blastout\\pssm',11)
   

    
                    
                