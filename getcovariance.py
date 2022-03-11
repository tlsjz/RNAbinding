import os
import linecache
import pickle
from runblast import psiblast

class getcovariance():
    def runMview(self,psiblastPath,mviewPath):
        psiblastFile = os.listdir(psiblastPath)
        for eachFile in psiblastFile:
            uniprotID = eachFile.split('.')[0]
            cmd='perl /home/songjiazhi/mview-1.64/bin/mview -in blast -out fasta '+psiblastPath+'/'+eachFile+' > '+mviewPath+'/'+uniprotID+'.aln.fasta'
            os.system(cmd)
    
    def MviewParser(self,mviewPath,mviewParsedPath):
        mviewFile = os.listdir(mviewPath)
        for eachFile in mviewFile:
            lines=linecache.getlines(mviewPath+'/'+eachFile)
            count = len(lines)
            i = 0
            parsedFile = open(mviewParsedPath+'/'+eachFile,'w')
            while i <= count-1:
                if i==count-1:
                    parsedFile.write(lines[i])
                else:
                    if lines[i][0]!='>' and lines[i+1][0]!='>':
                        newline=lines[i].replace('\n','')
                        lines[i]=newline
                    parsedFile.write(lines[i])
                i=i+1 
            parsedFile.close()
            
    def FormatCovarianceInput(self,mviewParsedPath,covarianceInputPath):
        parsedFile = os.listdir(mviewParsedPath)
        for eachFile in parsedFile:
            uniportID = eachFile.split('.')[0]
            lines = linecache.getlines(mviewParsedPath+'/'+eachFile)
            count = len(lines)
            i = 0
            inputFile = open(covarianceInputPath+'/'+uniportID+'.txt','w')
            while i < count-1:
                part1 = lines[i].split()[0][1:]
                part2 = lines[i+1]
                line = part1+'        '+part2
                inputFile.write(line)
                i = i+2
            inputFile.close()
    
    
    def runCovariance(self,method,covarianceInputPath,covarianceOutputPath):
        covarianceInputFile = os.listdir(covarianceInputPath)
        for eachFile in covarianceInputFile:
            uniprotID = eachFile.split('.')[0]
            cmd = 'java covariance.algorithms.'+method+'Covariance '+covarianceInputPath+'/'+eachFile+' '+covarianceOutputPath+'/'+eachFile
            os.system(cmd)
            
    def CovarianceParser(self,method,covarianceOutputPath):
        covarianceOutputFile = os.listdir(covarianceOutputPath)
        covarianceDic = {}
        for eachFile in covarianceOutputFile:
            uniprotID = eachFile.split('.')[0]
            #print(uniprotID)
            covarianceDic[uniprotID]={}
            lines = linecache.getlines(covarianceOutputPath+'/'+eachFile)
            for line in lines[1:]:
                residue1 = line.split()[0]
                residue2 = line.split()[1]
                covarianceValue = line.split()[2]
                covarianceDic[uniprotID][residue1+'_'+residue2] = float(covarianceValue)
        return covarianceDic
    
    def main(self,flag,method,fastaPath,outPath,pssmPath,mviewPath,mviewParsedPath,covarianceInputPath,covarianceOutputPath):
        if flag == '0':
            a = psiblast()
            a.runblast(fastaPath,outPath,pssmPath)
            self.runMview(outPath,mviewPath)
            self.MviewParser(mviewPath,mviewParsedPath)
            self.FormatCovarianceInput(mviewParsedPath,covarianceInputPath)
            self.runCovariance(method,covarianceInputPath,covarianceOutputPath)
            result = self.CovarianceParser(method,covarianceOutputPath)
            return result
        elif flag == '1':
            self.runMview(outPath,mviewPath)
            self.MviewParser(mviewPath,mviewParsedPath)
            self.FormatCovarianceInput(mviewParsedPath,covarianceInputPath)
            self.runCovariance(method,covarianceInputPath,covarianceOutputPath)
            result = self.CovarianceParser(method,covarianceOutputPath) 
            return result
        elif flag == '2':
            self.runCovariance(method,covarianceInputPath,covarianceOutputPath)
            result = self.CovarianceParser(method,covarianceOutputPath)
            return result
        elif flag == '3':
            result = self.CovarianceParser(method,covarianceOutputPath)
            return result
        elif flag == '4':
            self.runMview(outPath,mviewPath)
            self.MviewParser(mviewPath,mviewParsedPath)
            self.FormatCovarianceInput(mviewParsedPath,covarianceInputPath)
            self.runCovariance(method,covarianceInputPath,covarianceOutputPath)
        else:
            print('wrong parameter!')
            

if __name__=="__main__":
    co = getcovariance()
    CovarianceDic = co.main('1','ELSC','/home/songjiazhi/dnabinding/fasta','/home/songjiazhi/dnabinding/blastout/out','/home/songjiazhi/atpbinding/dnabinding/pssm','/home/songjiazhi/dnabinding/mviewout/original','/home/songjiazhi/dnabinding/mviewout/parsed','/home/songjiazhi/dnabinding/mviewout/covarianceinput','/home/songjiazhi/dnabinding/covarianceout')
    #print(a)
    pickleFile = open('/home/songjiazhi/dnabinding/ELSCdic.pickle','wb')
    pickle.dump(CovarianceDic,pickleFile)
            
        
                
                
        
        
    
    

                
        
            
            
            
        