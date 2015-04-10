__author__ = "Sindhu Ghanta"
__copyright__ = "Copyright 2015, Becklab"
__license__ = "Becklab"
__version__ = "1.0.0"
__maintainer__ = "Sindhu Ghanta"
__email__ = "sghanta2@bidmc.harvard.edu"
__status__ = "Completed"

''' This code reads all the *.jpg images from a given path and extracts brown 
    stain pixels in Epithelium and Stroma. Output is a binary image with 
    epi-stroma segmentation and a text file with corresponding number of pixels. 
    Input: Path to the folder
    Output: Binary image file
            Text file 
    Example Usage: brown_stain_extraction.py 'C:/Images/'
    This code works the best when image resolution is approximately 3 mu m/pixel.
'''

import sys
import glob
import os
from skimage import io, segmentation, util
import numpy as np
from skimage.feature import greycoprops, greycomatrix, local_binary_pattern
import pickle

def getNucleusPixel(RedArray,GreenArray,BlueArray):
    imageGreen_new = GreenArray[RedArray<float(150)/255]
    imageBlue_new = BlueArray[RedArray<float(150)/255]
    imageBlue_new1 = imageBlue_new[imageGreen_new<float(150)/255]
    imageBlue_new2 = imageBlue_new1[imageBlue_new1>float(150)/255]
    return len(imageBlue_new2)
    
def getBrownPixel(RedArray,GreenArray,BlueArray):
    imageGreen_new = GreenArray[RedArray>float(100)/255]
    imageBlue_new = BlueArray[RedArray>float(100)/255]
    imageBlue_new1 = imageBlue_new[imageGreen_new<float(100)/255]
    imageBlue_new2 = imageBlue_new1[imageBlue_new1<float(100)/255]
    return len(imageBlue_new2)
    
def getBackgroundPixel(RedArray,GreenArray,BlueArray):
    PIXEL_INTENSITY_SUM_THRESH = 0.9
    imageGreen_new = GreenArray[RedArray>PIXEL_INTENSITY_SUM_THRESH]
    imageBlue_new = BlueArray[RedArray>PIXEL_INTENSITY_SUM_THRESH]
    imageBlue_new1 = imageBlue_new[imageGreen_new>PIXEL_INTENSITY_SUM_THRESH]
    imageBlue_new2 = imageBlue_new1[imageBlue_new1>PIXEL_INTENSITY_SUM_THRESH]
    return len(imageBlue_new2) 
    
   
# This threshold is on the intensity value of each color channel beyong which it is considered background
PIXEL_INTENSITY_SUM_THRESH = 0.8
FRACTION_ONES_THRESH = 0.4
def processImage(imagePath,clf):
    
    # Read the image  
    color_original = util.img_as_float(io.imread(imagePath))+float(1)/255
    # Extract the r, g, b components of the image
    image_originalR, image_originalG, image_originalB = color_original[:,:,0], color_original[:,:,1], color_original[:,:,2]
    
              
    ###########################################################################
    # Initialize the texture components of the three channels in the image
    correlation_meanR = np.zeros((image_originalR.shape[0],image_originalR.shape[1]))
    correlation_meanG = np.zeros((image_originalG.shape[0],image_originalG.shape[1]))
    correlation_meanB = np.zeros((image_originalB.shape[0],image_originalB.shape[1]))
    correlation_varianceR = np.zeros((image_originalR.shape[0],image_originalR.shape[1]))
    correlation_varianceG = np.zeros((image_originalG.shape[0],image_originalG.shape[1]))
    correlation_varianceB = np.zeros((image_originalB.shape[0],image_originalB.shape[1]))
    contrast_meanR = np.zeros((image_originalR.shape[0],image_originalR.shape[1]))
    contrast_meanG = np.zeros((image_originalG.shape[0],image_originalG.shape[1]))
    contrast_meanB = np.zeros((image_originalB.shape[0],image_originalB.shape[1]))
    contrast_varianceR = np.zeros((image_originalR.shape[0],image_originalR.shape[1]))
    contrast_varianceG = np.zeros((image_originalG.shape[0],image_originalG.shape[1]))
    contrast_varianceB = np.zeros((image_originalB.shape[0],image_originalB.shape[1]))
    dissimilarity_meanR = np.zeros((image_originalR.shape[0],image_originalR.shape[1]))
    dissimilarity_meanG = np.zeros((image_originalG.shape[0],image_originalG.shape[1]))
    dissimilarity_meanB = np.zeros((image_originalB.shape[0],image_originalB.shape[1]))
    dissimilarity_varianceR = np.zeros((image_originalR.shape[0],image_originalR.shape[1]))
    dissimilarity_varianceG = np.zeros((image_originalG.shape[0],image_originalG.shape[1]))
    dissimilarity_varianceB = np.zeros((image_originalB.shape[0],image_originalB.shape[1]))
    homogeneity_meanR = np.zeros((image_originalR.shape[0],image_originalR.shape[1]))
    homogeneity_meanG = np.zeros((image_originalG.shape[0],image_originalG.shape[1]))
    homogeneity_meanB = np.zeros((image_originalB.shape[0],image_originalB.shape[1]))
    homogeneity_varianceR = np.zeros((image_originalR.shape[0],image_originalR.shape[1]))
    homogeneity_varianceG = np.zeros((image_originalG.shape[0],image_originalG.shape[1]))
    homogeneity_varianceB = np.zeros((image_originalB.shape[0],image_originalB.shape[1]))
    localBinaryPattern_meanR = np.zeros((image_originalR.shape[0],image_originalR.shape[1]))
    localBinaryPattern_meanG = np.zeros((image_originalG.shape[0],image_originalG.shape[1]))
    localBinaryPattern_meanB = np.zeros((image_originalB.shape[0],image_originalB.shape[1]))
    localBinaryPattern_varianceR = np.zeros((image_originalR.shape[0],image_originalR.shape[1]))
    localBinaryPattern_varianceG = np.zeros((image_originalG.shape[0],image_originalG.shape[1]))
    localBinaryPattern_varianceB = np.zeros((image_originalB.shape[0],image_originalB.shape[1]))
    
    # Calculate the texture components
    for row in range(19,image_originalR.shape[0]):               
        for col in range(19,image_originalR.shape[1]):            
            image_patchR = image_originalR[row-19:row,col-19:col]
            image_patchG = image_originalG[row-19:row,col-19:col]
            image_patchB = image_originalB[row-19:row,col-19:col]
            
            gR = greycomatrix(np.uint8(image_patchR*49), [1,5], [0, np.pi/2], 
                                levels=50,normed=True, symmetric=True)
            gG = greycomatrix(np.uint8(image_patchG*49), [1,5], [0, np.pi/2], 
                                levels=50,normed=True, symmetric=True)
            gB = greycomatrix(np.uint8(image_patchB*49), [1,5], [0, np.pi/2], 
                                levels=50,normed=True, symmetric=True)                
            correlationR = greycoprops(gR, 'correlation')                
            correlation_meanR[row-9,col-9] = correlationR.mean()
            correlation_varianceR[row-9,col-9] = correlationR.var()
            correlationG = greycoprops(gG, 'correlation')
            correlation_meanG[row-9,col-9] = correlationG.mean()
            correlation_varianceG[row-9,col-9] = correlationG.var()
            correlationB = greycoprops(gB, 'correlation')
            correlation_meanB[row-9,col-9] = correlationB.mean()
            correlation_varianceB[row-9,col-9] = correlationB.var()
                
            contrastR = greycoprops(gR,'contrast')
            contrastG = greycoprops(gG,'contrast')
            contrastB = greycoprops(gB,'contrast')
            contrast_meanR[row,col] = contrastR.mean()
            contrast_varianceR[row-9,col-9] = contrastR.var()
            contrast_meanG[row-9,col-9] = contrastG.mean()
            contrast_varianceG[row-9,col-9] = contrastG.var()
            contrast_meanB[row-9,col-9] = contrastB.mean()
            contrast_varianceB[row-9,col-9] = contrastB.var()
                
            dissimilarityR = greycoprops(gR,'dissimilarity')
            dissimilarityG = greycoprops(gG,'dissimilarity')
            dissimilarityB = greycoprops(gB,'dissimilarity')
            dissimilarity_meanR[row-9,col-9] = dissimilarityR.mean()
            dissimilarity_varianceR[row-9,col-9] = dissimilarityR.var() 
            dissimilarity_meanG[row-9,col-9] = dissimilarityG.mean()
            dissimilarity_varianceG[row-9,col-9] = dissimilarityG.var() 
            dissimilarity_meanB[row-9,col-9] = dissimilarityB.mean()
            dissimilarity_varianceB[row-9,col-9] = dissimilarityB.var() 
                
            homogeneityR = greycoprops(gR,'homogeneity')
            homogeneityG = greycoprops(gG,'homogeneity')
            homogeneityB = greycoprops(gB,'homogeneity')
            homogeneity_meanR[row-9,col-9] = homogeneityR.mean()
            homogeneity_varianceR[row-9,col-9] = homogeneityR.var()  
            homogeneity_meanG[row-9,col-9] = homogeneityG.mean()
            homogeneity_varianceG[row-9,col-9] = homogeneityG.var() 
            homogeneity_meanB[row-9,col-9] = homogeneityB.mean()
            homogeneity_varianceB[row-9,col-9] = homogeneityB.var() 
                
            localBinaryPatternR = local_binary_pattern(image_patchR,8,3)  
            localBinaryPatternG = local_binary_pattern(image_patchG,8,3) 
            localBinaryPatternB = local_binary_pattern(image_patchB,8,3) 
            localBinaryPattern_meanR[row-9,col-9] = localBinaryPatternR.mean()
            localBinaryPattern_varianceR[row-9,col-9] = localBinaryPatternR.var()
            localBinaryPattern_meanG[row-9,col-9] = localBinaryPatternG.mean()
            localBinaryPattern_varianceG[row-9,col-9] = localBinaryPatternG.var()
            localBinaryPattern_meanB[row-9,col-9] = localBinaryPatternB.mean()
            localBinaryPattern_varianceB[row-9,col-9] = localBinaryPatternB.var()
            
       
    # Normalize the calculated texture components to [0,1]
    correlation_meanR = correlation_meanR/correlation_meanR.max()
    correlation_varianceR = correlation_varianceR/correlation_varianceR.max()
    correlation_meanG = correlation_meanG/correlation_meanG.max()
    correlation_varianceG = correlation_varianceG/correlation_varianceG.max()
    correlation_meanB = correlation_meanB/correlation_meanB.max()
    correlation_varianceB = correlation_varianceB/correlation_varianceB.max()
    
    contrast_meanR = contrast_meanR/contrast_meanR.max()
    contrast_varianceR = contrast_varianceR/contrast_varianceR.max()
    contrast_meanG = contrast_meanG/contrast_meanG.max()
    contrast_varianceG = contrast_varianceG/contrast_varianceG.max()
    contrast_meanB = contrast_meanB/contrast_meanB.max()
    contrast_varianceB = contrast_varianceB/contrast_varianceB.max()
    
    dissimilarity_meanR = dissimilarity_meanR/dissimilarity_meanR.max()
    dissimilarity_varianceR = dissimilarity_varianceR/dissimilarity_varianceR.max()
    dissimilarity_meanG = dissimilarity_meanG/dissimilarity_meanG.max()
    dissimilarity_varianceG = dissimilarity_varianceG/dissimilarity_varianceG.max()
    dissimilarity_meanB = dissimilarity_meanB/dissimilarity_meanB.max()
    dissimilarity_varianceB = dissimilarity_varianceB/dissimilarity_varianceB.max()
    
    homogeneity_meanR = homogeneity_meanR/homogeneity_meanR.max()
    homogeneity_varianceR = homogeneity_varianceR/homogeneity_varianceR.max()
    homogeneity_meanG = homogeneity_meanG/homogeneity_meanG.max()
    homogeneity_varianceG = homogeneity_varianceG/homogeneity_varianceG.max()
    homogeneity_meanB = homogeneity_meanB/homogeneity_meanB.max()
    homogeneity_varianceB = homogeneity_varianceB/homogeneity_varianceB.max()
    
    localBinaryPattern_meanR = localBinaryPattern_meanR/localBinaryPattern_meanR.max()
    localBinaryPattern_varianceR = localBinaryPattern_varianceR/localBinaryPattern_varianceR.max()
    localBinaryPattern_meanG = localBinaryPattern_meanG/localBinaryPattern_meanG.max()
    localBinaryPattern_varianceG = localBinaryPattern_varianceG/localBinaryPattern_varianceG.max()
    localBinaryPattern_meanB = localBinaryPattern_meanB/localBinaryPattern_meanB.max()
    localBinaryPattern_varianceB = localBinaryPattern_varianceB/localBinaryPattern_varianceB.max()        
    
    
    ############################################################################
    # Divide the image into a number of superpixels   
    numSegments = 800
    segments = segmentation.slic(color_original, n_segments = numSegments, compactness=6, max_iter=50, 
                                sigma = 5, enforce_connectivity=True)
    imageShape = color_original.shape  
     
    # Create a color image which is all black, let stroma be Red and Epithelium be Green
    SegStromaEpithelium = np.zeros((imageShape[0], imageShape[1], imageShape[2]))   
    
    # Initialize the statistics of stain
    A_Brown_St =0     
    A_Brown_Ep =0   
    A_St = 0   
    A_Ep = 0 
    N_St = 0
    N_Ep = 0
    
    # Calculate the class of each superpixel
    for segmentNum in range(0,segments.max()):
        newArray = color_original[segments==segmentNum] 
        
        # Extract the texture features of this superpixel
        corrMeanR = correlation_meanR[segments==segmentNum].mean()
        corrMeanG = correlation_meanG[segments==segmentNum].mean()
        corrMeanB = correlation_meanB[segments==segmentNum].mean()
        corrVarR =  correlation_varianceR[segments==segmentNum].mean()  
        corrVarG =  correlation_varianceG[segments==segmentNum].mean()
        corrVarB =  correlation_varianceB[segments==segmentNum].mean()
        
        contMeanR = contrast_meanR[segments==segmentNum].mean()
        contMeanG = contrast_meanG[segments==segmentNum].mean()
        contMeanB = contrast_meanB[segments==segmentNum].mean()
        contVarR = contrast_varianceR[segments==segmentNum].mean()
        contVarG = contrast_varianceG[segments==segmentNum].mean()
        contVarB = contrast_varianceB[segments==segmentNum].mean()
        
        dissMeanR = dissimilarity_meanR[segments==segmentNum].mean()
        dissMeanG = dissimilarity_meanG[segments==segmentNum].mean()
        dissMeanB = dissimilarity_meanB[segments==segmentNum].mean()
        dissVarR = dissimilarity_varianceR[segments==segmentNum].mean()
        dissVarG = dissimilarity_varianceG[segments==segmentNum].mean()
        dissVarB = dissimilarity_varianceB[segments==segmentNum].mean()
        
        homoMeanR = homogeneity_meanR[segments==segmentNum].mean()
        homoMeanG = homogeneity_meanG[segments==segmentNum].mean()
        homoMeanB = homogeneity_meanB[segments==segmentNum].mean()
        homoVarR = homogeneity_varianceR[segments==segmentNum].mean()
        homoVarG = homogeneity_varianceG[segments==segmentNum].mean()
        homoVarB = homogeneity_varianceB[segments==segmentNum].mean()
        
        lbpMeanR = localBinaryPattern_meanR[segments==segmentNum].mean()
        lbpMeanG = localBinaryPattern_meanG[segments==segmentNum].mean()
        lbpMeanB = localBinaryPattern_meanB[segments==segmentNum].mean()
        
        lbpVarR = localBinaryPattern_varianceR[segments==segmentNum].mean()
        lbpVarG = localBinaryPattern_varianceG[segments==segmentNum].mean()
        lbpVarB = localBinaryPattern_varianceB[segments==segmentNum].mean()
        
        RedArray = np.trim_zeros(newArray[:,0])
        GreenArray = np.trim_zeros(newArray[:,1])
        BlueArray = np.trim_zeros(newArray[:,2])         
        fractionOnes = float(getBackgroundPixel(RedArray,GreenArray,BlueArray))/len(RedArray)        
        indexPixels = np.where(segments==segmentNum)
        
        # If fraction of white pixels is too high, assume that it is background or fat and leave it out without processing
        if(fractionOnes<FRACTION_ONES_THRESH):
            classLabel = clf.predict([corrMeanR,corrMeanG,corrMeanB,corrVarR,corrVarG,corrVarB,
            contMeanR,contMeanG,contMeanB,contVarR,contVarG,contVarB,
            dissMeanR,dissMeanG,dissMeanB,dissVarR,dissVarG,dissVarB,
            homoMeanR,homoMeanG,homoMeanB,homoVarR,homoVarG,homoVarB,
            lbpMeanR,lbpMeanG,lbpMeanB,lbpVarR,lbpVarG,lbpVarB])    
            if(classLabel == 1): # Its epithelium
                SegStromaEpithelium[indexPixels[0],indexPixels[1],0] = 1 
                A_Brown_Ep = A_Brown_Ep + getBrownPixel(RedArray,GreenArray,BlueArray)
                A_Ep = A_Ep + len(RedArray)  
                N_Ep = N_Ep + getNucleusPixel(RedArray,GreenArray,BlueArray)
            else:  # Its Stroma
                SegStromaEpithelium[indexPixels[0],indexPixels[1],1] = 1 
                A_Brown_St = A_Brown_St + getBrownPixel(RedArray,GreenArray,BlueArray)
                A_St = A_St + len(RedArray)  
                N_St = N_St + getNucleusPixel(RedArray,GreenArray,BlueArray)
    
    return (SegStromaEpithelium,A_St,A_Ep,A_Brown_St,A_Brown_Ep,N_St,N_Ep)



def processFolder(FolderPath):
    path = FolderPath  
     
    # Create a text file to save the results
    f = open(path+"Results.txt","w")
    f.write("Image path , Stroma pixels, Epithelium pixels, Brown Spots in Stroma, Brown Spots in Epithelium, Nucleus Pixels in Stroma, Nucleus pixels in Epithelium \n")
    
    # Change the path to directory containing images
    os.chdir(path)
    
    # Load the SVM parameters
    clf = pickle.load(open('RGB_SVM.pk1'))
    
    # Process all the JPEG images in the folder given by user   
    for file in glob.glob("*.jpg"):
        imagePath = path+file
        
        # Call the function to process the image
        [SegStromaEpithelium,A_St,A_Ep,A_Brown_St,A_Brown_Ep,N_St,N_Ep] = processImage(imagePath,clf)
        print("Processing " + file)
        
        # Save the segmented image
        io.imsave(path+file[0:len(file)-4]+"_binarySegementation.jpg",SegStromaEpithelium)
        
        # Save the brown stain statistics
        f.write(imagePath+ " , "+ str(A_St)+ " , " + str(A_Ep) + " , " + str(A_Brown_St)+" , "+str(A_Brown_Ep)+ " , " + str(N_St)+" , "+str(N_Ep)+"\n")
    
    # Close the text file after processing is finished
    f.close()

        
if __name__ == "__main__":
   processFolder(sys.argv[1])

        
