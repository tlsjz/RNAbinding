import os
import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
#import matplotlib.pyplot as plt

modelpickle = open('D:\\RNAbinding\\RB344\\newfeature\\model\\modeleleven.model','rb')
model = pickle.load(modelpickle)
filelist = os.listdir('D:\\RNAbinding\\RB75\\pdbfasta')
testlabel = []
predictlabel = []
template = 0
featurenum = 0
overlap = 0
mccdic = {}
flist = []
#site3pickle = open('D:\\RNAbinding\\RB75\\3.5\\sitedic.pickle','rb')
#site3 = pickle.load(site3pickle)
for eachfile in filelist:
    fullid = eachfile.split('.')[0]
    #print(fullid)
    pdbid = fullid.split('_')[0]
    chainid = fullid.split('_')[1]
    inputFilePickle = open('D:\\RNAbinding\\RB75\\inputfile\\'+fullid+'.pickle','rb')
    inputfile = pickle.load(inputFilePickle)
    eachlabel = inputfile[0]
    feature = inputfile[1]
    Fscorelist = []
    FeaturePredict = model.predict_proba(feature)
    for each in FeaturePredict:
        Fscorelist.append(each[1])
    length = len(eachlabel)
    eachpredict = []      
    for i in range(0,length):
        #eachpredict.append(Fscorelist[i])
        if Fscorelist[i] > 0.5:
            eachpredict.append(Fscorelist[i])
        else:
            if i < 4:
                proba = Fscorelist[i]
                if Fscorelist[i+1] > 0.5:
                    proba = proba+Fscorelist[i+1]*(1/2)
                elif Fscorelist[i+2] > 0.5:
                    proba = proba+Fscorelist[i+2]*(1/4)
                elif Fscorelist[i+3] > 0.5:
                    proba = proba+Fscorelist[i+3]*(1/8)
                elif Fscorelist[i+4] > 0.5:
                    proba = proba+Fscorelist[i+4]*(1/16)
                #elif Fscorelist[i+5] > 0.5:
                    #proba = proba+Fscorelist[i+5]*(1/32)
                #elif Fscorelist[i+6] > 0.5:
                    #proba = proba+Fscorelist[i+6]*(1/64)
                eachpredict.append(proba)
            elif i > len(eachlabel)-5:
                proba = Fscorelist[i]
                if Fscorelist[i-1] > 0.5:
                    proba = proba+Fscorelist[i-1]*(1/2)
                elif Fscorelist[i-2] > 0.5:
                    proba = proba+Fscorelist[i-2]*(1/4)
                elif Fscorelist[i-3] > 0.5:
                    proba = proba+Fscorelist[i-3]*(1/8)
                elif Fscorelist[i-4] > 0.5:
                    proba = proba+Fscorelist[i-4]*(1/16)
                #elif Fscorelist[i-5] > 0.5:
                    #proba = proba+Fscorelist[i-5]*(1/32)
                #elif Fscorelist[i-6] > 0.5:
                    #proba = proba+Fscorelist[i-6]*(1/64)
                eachpredict.append(proba)      
            else:
                proba = Fscorelist[i]
                if Fscorelist[i-1] > 0.5 or Fscorelist[i+1] > 0.5:
                    if Fscorelist[i-1] >= Fscorelist[i+1]:
                        proba = proba+Fscorelist[i-1]*(1/2)
                    else:
                        proba = proba+Fscorelist[i+1]*(1/2)
                elif Fscorelist[i+2] > 0.5 or Fscorelist[i-2] > 0.5:
                    if Fscorelist[i-2] >= Fscorelist[i+2]:
                        proba = proba+Fscorelist[i-2]*(1/4)
                    else:
                        proba = proba+Fscorelist[i+2]*(1/4)
                elif Fscorelist[i+3] > 0.5 or Fscorelist[i-3] > 0.5:
                    if Fscorelist[i-3] >= Fscorelist[i+3]:
                        proba = proba+Fscorelist[i-3]*(1/8)
                    else:
                        proba = proba+Fscorelist[i+3]*(1/8)  
                elif Fscorelist[i+4] > 0.5 or Fscorelist[i-4] > 0.5:
                    if Fscorelist[i-4] >= Fscorelist[i+4]:
                        proba = proba+Fscorelist[i-4]*(1/16)
                    else:
                        proba = proba+Fscorelist[i+4]*(1/16) 
                #elif Fscorelist[i+5] > 0.5 or Fscorelist[i-5] > 0.5:
                    #if Fscorelist[i-5] >= Fscorelist[i+5]:
                        #proba = proba+Fscorelist[i-5]*(1/32)
                    #else:
                        #proba = proba+Fscorelist[i+5]*(1/32)
                #elif Fscorelist[i+6] > 0.5 or Fscorelist[i-6] > 0.5:
                    #if Fscorelist[i-6] >= Fscorelist[i+6]:
                        #proba = proba+Fscorelist[i-6]*(1/64)
                    #else:
                        #proba = proba+Fscorelist[i+6]*(1/64)                 
                eachpredict.append(proba)  
    eachresult = []
    for i in range(0,length):
        #eachresult.append(eachpredict[i])
        if eachpredict[i] > 0.5:
            eachresult.append(1)
        else:
            eachresult.append(0)
            
    templatePickle = open('D:\\RNAbinding\\RB75\\templatePrediction\\'+fullid+'.pickle','rb')
    templatePredict = pickle.load(templatePickle)
    Tscorelist = []
    if templatePredict['score'] > 50:
        for i in range(0,length):
            if i in templatePredict['SeqId']:
                Tscorelist.append(1)
            else:
                Tscorelist.append(0)
    else:
        for i in range(0,length):
            Tscorelist.append(0)
            
    for i in range(0,length):
        if eachpredict[i] >= 0.5 or Tscorelist[i] == 1:
        #if Tscorelist[i] == 1:
            eachresult.append(1)
        else:
            eachresult.append(0)
            
    #if fullid == '3izv_M':
        #print(eachlabel)
        #print(eachresult)
        #print(metrics.matthews_corrcoef(eachlabel,eachresult))
        
    #for i in range(0,length):
        #if eachpredict[i] >= 0.5 and Tscorelist[i] == 1:
            #overlap = overlap+1
        #elif eachpredict[i] >= 0.5 and Tscorelist[i] != 1:
            #featurenum = featurenum+1
        #elif Tscorelist[i] == 1 and eachpredict[i] < 0.5:
            #template = template+1
    #print(metrics.matthews_corrcoef(eachlabel,eachresult))    
    #mccdic[fullid] = metrics.matthews_corrcoef(eachlabel,eachresult)
    #flist.append(metrics.f1_score(eachlabel,eachresult))
    for each in eachlabel:
        testlabel.append(each)
    for each in eachresult:
        predictlabel.append(each)
    #for i in range(length):
        #print(eachlabel[i],Fscorelist[i],eachpredict[i])
    #for i in range(0,length):
        #testlabel.append(eachlabel[i])
    #for i in range(0,length):
        #predictlabel.append(Fscorelist[i])
#plt.figure()
#lw = 2
print(testlabel)
print(predictlabel)
#print(mccdic)
#fpr,tpr,threshold = metrics.roc_curve(testlabel,predictlabel)
#roc_auc = metrics.auc(fpr,tpr)
#plt.plot(fpr,tpr,color = 'darkorange',lw = lw,label = 'RF ROC curve(area = %0.2f)' %roc_auc)
#plt.xlim([0.0,1.0])
#plt.ylim([0.0,1.0])
#plt.grid()
#plt.xlabel('false positive label')
#plt.ylabel('true positive label')
#plt.title('ROC')
#plt.legend(loc = 'lower right')
#plt.show()
#print(roc_auc)
print(metrics.accuracy_score(testlabel,predictlabel))
print(metrics.precision_score(testlabel,predictlabel))
print(metrics.recall_score(testlabel,predictlabel))
print(metrics.matthews_corrcoef(testlabel,predictlabel))
print(metrics.f1_score(testlabel,predictlabel))
print(overlap,template,featurenum)
#mcctotal = 0
#for each in mcclist:
    #mcctotal = mcctotal+each
#print(mcctotal/len(mcclist))
#ftotal = 0
#for each in flist:
    #ftotal = ftotal+each
#print(ftotal/len(flist))

            
    
        
    