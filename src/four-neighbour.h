#ifndef FOURNEIGHBOUR_H
#define FOURNEIGHBOUR_H
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include <vector>
#include <utility>
#include <algorithm>
#include <iostream>

typedef unsigned char ScalarPixelType;
typedef itk::Image <ScalarPixelType, 2> ScalarImageType;
typedef itk::ImageFileReader<ScalarImageType> ScalarReaderType;
typedef itk::ImageRegionIterator<ScalarImageType> IteratorType;
void recordRuns(ScalarImageType *input, 
                std::vector<int> &cstRun, std::vector<int> &cenRun, 
                std::vector<int> &rRun){
    IteratorType iterator(input, input->GetRequestedRegion());
    iterator.GoToBegin();
    int row,column = 1;
    int rsize = 512, csize = 512;
    ScalarPixelType curPixel = iterator.Get();
    ScalarPixelType prePixel;
    for(row = 1;row<=rsize;row++){
        for(column = 1;column <=csize;iterator++,column++){
            prePixel = curPixel;
            curPixel = iterator.Get();
            if(column ==1&&iterator.Get()==0){
                rRun.push_back(row);
                cstRun.push_back(1);
            }
            else if(column>1){
                if(prePixel >0 &&curPixel ==0){
                    rRun.push_back(row);
                    cstRun.push_back(column);
                }
                if(prePixel==0&&curPixel >0){
                    //rRun.push_back(row);
                    cenRun.push_back(column-1);
                }
                if(column == csize &&curPixel==0){
                    //rRun.push_back(row);
                    cenRun.push_back(column);
                }
            }
        }
    }
}

void firstPass(std::vector<int> &cstRun, std::vector<int> &cenRun, std::vector<int> &rRun, std::vector<int> &tag,
               std::vector<std::pair<int,int>> &equivalence){
    int tagIndex = 1;
    int curRow = 1;
    std::vector<int> preSt, preEn, preTag, curSt, curEn, curTag;
    for(int i=0;i<cstRun.size();i++){
        if(rRun[i] == 1){
            tag.push_back(tagIndex);
            curSt.push_back(cstRun[i]);
            curEn.push_back(cenRun[i]);
            curTag.push_back(tagIndex);
            tagIndex++;
        }
        else{
            if(rRun[i] != curRow){
                curRow = rRun[i];
                preSt = curSt;
                preEn = curEn;
                preTag = curTag;
                curSt.clear();
                curEn.clear();
                curTag.clear();
            }
            curSt.push_back(cstRun[i]);
            curEn.push_back(cenRun[i]);
            std::vector<int> intersectTag;
            for(int j=0;j<preSt.size();j++){
                if(preEn[j]>=cstRun[i] && preSt[j]<= cenRun[i]){
                    intersectTag.push_back(preTag[j]);
                }
            }
            //intersectTag是否为空
            if(!intersectTag.size()){
                tag.push_back(tagIndex);
                curTag.push_back(tagIndex);
                tagIndex++;
            }
            else{
                //取小并记录等价对
                //tag.push_back(*std::min_element(std::begin(intersectTag), std::end(intersectTag)));
                std::sort(std::begin(intersectTag), std::end(intersectTag));
                tag.push_back(intersectTag[0]);
                curTag.push_back(intersectTag[0]);
                for(int k=1;k<intersectTag.size();k++){
                    equivalence.push_back(std::make_pair(intersectTag[0],intersectTag[k]));
                }
            }
        }
    }
}

//合并等价对
void secondPass(std::vector<int> &tag, std::vector<std::pair<int,int>> &equivalence){
    int maxTag = *std::max_element(tag.begin(),tag.end());
    //std::cout<<"maxtag:"<<maxTag<<std::endl;
    std::vector<int> remainTag,tempList;
    std::vector<bool> tagFlag(maxTag,false);
    std::vector<std::pair<int,int>> replaceTable;
    bool flag = false;
    int maxDepth = maxTag;
    for(int i=maxTag;i>=1;i--)
        remainTag.push_back(i);
    for(std::vector<int>::iterator iter = remainTag.begin();iter!=remainTag.end();iter++)
    {
        flag = false;
        maxDepth = *iter;
        while(!flag){
            flag = true;
            //tempList.push_back(*iter);
            for(int j = 0;j<equivalence.size();j++){
                if(equivalence[j].second == maxDepth){
                    if(equivalence[j].first<maxDepth){
                        maxDepth = equivalence[j].first;
                        flag = false;
                    }
                    //tempList.push_back(equivalence[j].first);
                }
            }
        }
        if(flag){
            replaceTable.push_back(std::make_pair(*iter,maxDepth));
        }
    }
    //std::cout<<"replaceTable:"<<replaceTable.size()<<std::endl;
    //std::cout<<"maxTag:"<<maxTag<<std::endl;
    //replaceTable.size() == maxTag;
    for(int i=0;i<tag.size();i++){
        tag[i] = replaceTable[maxTag-tag[i]].second;
    }
}

void calculateArea(std::vector<int>&cstRun, std::vector<int>&cenRun, std::vector<int>&tag, std::vector<std::pair<int,int>> &sortedArea){
    //std::vector<bool> initialized(*std::max_element(tag.begin(),tag.end()),false);
    std::map<int,int> area;
    for(std::vector<int>::iterator iter = tag.begin();iter!=tag.end();iter++){
        area[*iter] = 0;
    }
    for(int i=0;i<tag.size();i++){
        area[tag[i]] += cenRun[i] - cstRun[i] +1;
    }
    std::vector<std::pair<int,int>> sorted(area.begin(),area.end());
    struct CmpByValue{
        bool operator()(const std::pair<int,int>&lhs,const std::pair<int,int>&rhs){
            return lhs.second>rhs.second;
        }
    };
    std::sort(sorted.begin(),sorted.end(),CmpByValue());
    sortedArea = sorted;
}
#endif