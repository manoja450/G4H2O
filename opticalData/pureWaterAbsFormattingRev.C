
void pureWaterAbsFormattingRev(){

    //will read in lambda vs. abs coeff (nm, cm-1 for this file)
    //will write out lambda (nm) vs. abs length (cm)
    
    TString inFile("opticalData/pureWaterAbsCoeff.dat");
    TString outFile("opticalData/pureWaterAbsLength_formattedrev.dat");
    
    Double_t minLambda = 200;
    Double_t maxLambda = 800;
    
    Int_t numEn = 251;
    TGraph *g = new TGraph(251);
    
    //https://omlc.org/spectra/water/data/buiteveld94.dat
    ifstream in(inFile);
    
    Double_t minL = 9999999, maxL = 0;
    
    //nm, cm-1
    Double_t lambda, coeff;
    for(Int_t iEn=0; iEn<numEn; iEn++){
        in >> lambda >> coeff;
        
        g->SetPoint(iEn,lambda,coeff);
        if(lambda>maxL) maxL = lambda;
        if(lambda<minL) minL = lambda;
    }
    in.close();
    
    Int_t numEn2 = numEn;

    if(minL > minLambda) numEn2++;
    if(maxL < maxLambda) numEn2++;

    //now create data with extrapolated values to
    //get to min and max lambda if necessary
    
    TGraph *g2 = new TGraph(numEn2);
    
    Int_t iEnTemp = 0;
    if(minL > minLambda){
        Double_t newVal = g->Eval(minLambda);
        if(newVal>0) newVal = 1.0/newVal;
        g2->SetPoint(iEnTemp++,minLambda,newVal);
    }

    Double_t xt, yt;
    for(Int_t iEn=0; iEn<numEn; iEn++){
        g->GetPoint(iEn,xt,yt);
        Double_t newVal = yt;
        if(newVal>0) newVal = 1.0/newVal;
        g2->SetPoint(iEnTemp++,xt,newVal);
        g->SetPoint(iEn,xt,newVal);
    }
    if(maxL < maxLambda){
        Double_t newVal = g->Eval(maxLambda);
        if(newVal>0) newVal = 1.0/newVal;
        g2->SetPoint(iEnTemp++,maxLambda,newVal);
    }
    
    g2->Draw("alp");
    g->SetLineColor(kRed);
    g->Draw("lp");
    
    FILE *fout = fopen(outFile,"w");
    fprintf(fout,"%d\n",numEn2);
    for(Int_t iEn=0; iEn<numEn2; iEn++){
        g2->GetPoint(iEn,xt,yt);
        fprintf(fout,"%.1f\t%e\n",xt,yt);
    }
    fclose(fout);

}
