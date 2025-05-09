#include <iostream>
#include <fstream>
#include <TGraphPolar.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TGraphPolargram.h>
#include <TMath.h>
#include <TLatex.h>
#include <TStyle.h>

std::string LabelsPolarTextual24Divs[24] = {
    "N",
    "15",
    "30",
    "NE",
    "60",
    "75",
    "E",
    "105",
    "120",
    "SE",
    "150",
    "165",
    "S",
    "195",
    "210",
    "SW",
    "240",
    "255",
    "W",
    "285",
    "300",
    "NW",
    "330",
    "345"
};

int main(){
	//Useful variables and conversions
	auto pi = TMath::Pi();
    double degToRad = pi / 180;
    
    //NOTE! The ROOT standard for polar graphing has north to the right.
    //Thus, the Azimuth 0 is parallel to positive X in ROOT's internal system.
    //Thanks, ROOT.
    
    double CelestialPoleY = 0; //  by definition
    double CelestialPoleX = 44.5892595; //Latitude of the observatory (thanks, google maps) :)
    
    //Root application initialization
    TApplication app("app", 0,0);   //initialization
    gStyle->SetOptTitle(0);			//removing title
    TCanvas c1;						//the canvas for graphing
    
    //Graph style choices
    auto polarFrame = new TGraphPolargram("polarFrame", 0, 360*degToRad, 0, 90);
    polarFrame->SetToRadian();
    polarFrame->SetNdivRadial(9);
    polarFrame->SetNdivPolar(24); //STANDARD IS 24. ANY OTHER DIV NUMBER WILL NOT HAVE LETTER LABELING
    polarFrame->Draw("P SAME");
    polarFrame->SetPolarLabelSize(0);
    double r = 1.1; // Relative distance from graph center, used to label the divisions.
    double directionOfNorth = 90; //Root's default is right.
    int PolarDivs = polarFrame->GetNdivPolar();
    int RadialDivs = polarFrame->GetNdivRadial();
    double divAngle = 360.0/PolarDivs; // Angle of a single polar division
    
    //Labelling of polar divisions.
    if(PolarDivs == 24){
	    for(int i = 0; i<PolarDivs; i++){
    	    double x = r * TMath::Cos((directionOfNorth + i * divAngle) * degToRad);
    	    double y = r * TMath::Sin((directionOfNorth + i * divAngle) * degToRad);
    	    TLatex *Label = new TLatex(x, y, LabelsPolarTextual24Divs[i].c_str());
    	    Label->SetTextAlign(22); // Center align
    	    Label->SetTextSize(0.05);
    	    Label->Draw();
    	}
    }
    else{
	    for(int i = 0; i<PolarDivs; i++){
    	    double x = r * TMath::Cos((directionOfNorth + i * divAngle) * degToRad);
    	    double y = r * TMath::Sin((directionOfNorth + i * divAngle) * degToRad);
    	    TLatex *Label = new TLatex(x, y, std::to_string(i).c_str());
    	    Label->SetTextAlign(22); // Center align
    	    Label->SetTextSize(0.05);
    	    Label->Draw();
    	}
    }
    
    //Labelling of radial division. 
    //In the telescope system, 90 is the horizon. In stellarium, 0 is the orizon. TODO: Display both with corresponding labelling 
    for(int i = 0; i<RadialDivs; i++){
        double y = (i+1) * 0.1 + 0.0125*i; //TODO: need better alignment, the division is non linear.
        double x = 0;
        TLatex *Label = new TLatex(x, y, std::to_string((i+1)*10).c_str());
        Label->SetTextAlign(25);
        Label->SetTextSize(0.02);
        Label->Draw();
    }
    
    //File handling
    
    std::ifstream takenFile("taken.dat");
    std::ifstream predictedFile("predicted.dat");
    
    double takenAz = 0, takenZd = 0;
    
    std::vector<double> takenAzList(0);
    std::vector<double> takenZdList(0);
    
    //Taken stars processing start
    
    while(takenFile >> takenAz >> takenZd){
    	takenAzList.push_back((takenAz + directionOfNorth)*degToRad); //Thanks, root. You really had to default north right.
    	takenZdList.push_back(takenZd);
    }
    
    //Brief explanation of the following few lines: ROOT is strange in its way of handling graphs (albeit really useful)
    //The _t type of variables is ROOT internal stuff, and its constructors work way better when using that type of variables.
    //Whenever a conversion is made that seems strange or useless, it's done to work around constructor limitations
    
    Double_t takenAzGraph[takenAzList.size()];  //Thanks, root
    Double_t takenZdGraph[takenZdList.size()];  //Thanks, root
    
    for(unsigned long int i=0; i<takenAzList.size(); i++){
        takenAzGraph[i] = takenAzList[i];
        takenZdGraph[i] = takenZdList[i];
    }
    
    //Graph creation and styling
    TGraphPolar gr1(takenAzList.size(), takenAzGraph, takenZdGraph);
    
    gr1.SetPolargram(polarFrame);
    gr1.SetMarkerStyle(29);
    gr1.SetMarkerColor(kBlack);
    gr1.SetMinRadial(0);
    gr1.SetMaxRadial(90);
    gr1.SetMarkerSize(1.5);
    gr1.Draw("P SAME");
    
    //Taken stars processing end
    
    //Prediction processing
    //This part of the tool is meant to be used to decide if an object will fall in its path inside a mapped sky zone.
    //It takes a pair of points as input, and calculates the path of the object from the starting point
    //The path calculation is independent from the end point to allow the user to give a single point (albeit repeated twice) to calculate.
    //TODO: change the input system to STARTPOINT+TIME instead of STARTPOINT+ENDPOINT
    //TODO: the formula is not mathematically accurate, the 
    
    double predictionAzStart = 0, predictionZdStart = 0, predictionAzEnd = 0, predictionZdEnd = 0;
    
    std::vector<double> predictionStartPointsAzList(0);
    std::vector<double> predictionStartPointsZdList(0);
    std::vector<double> predictionAzList(0);
    std::vector<double> predictionZdList(0);
    
    while(predictedFile >> predictionAzStart >> predictionZdStart >> predictionAzEnd >> predictionZdEnd ){
    
    	//The processing assumes that the user is looking for the coordinates on stellarium
    	//In stellarium, the horizon is 0
    	predictionStartPointsAzList.push_back((predictionAzStart) *degToRad);
    	predictionAzList.push_back((predictionAzStart + directionOfNorth) *degToRad);
    	predictionAzList.push_back((predictionAzEnd + directionOfNorth) *degToRad);
    	
    	predictionStartPointsZdList.push_back(90 - predictionZdStart);
    	predictionZdList.push_back(90 - predictionZdStart);
    	predictionZdList.push_back(90 - predictionZdEnd);
    }
    
    Double_t AzPredictedGraph[predictionAzList.size()];  //Thanks, root
    Double_t ZdPredictedGraph[predictionAzList.size()];  //Thanks, root


    for(unsigned long int i=0; i<predictionAzList.size(); i++){
        AzPredictedGraph[i] = predictionAzList[i];
        ZdPredictedGraph[i] = predictionZdList[i];
    }
    
    //Graph of the start and end points
    TGraphPolar gr2(predictionAzList.size(), AzPredictedGraph, ZdPredictedGraph);

    gr2.SetPolargram(polarFrame);
    gr2.SetMarkerStyle(29);
    gr2.SetMarkerColor(kRed);
    gr2.SetMinRadial(0);
    gr2.SetMaxRadial(90);
    gr2.SetMarkerSize(1.5);
    gr2.Draw("P SAME");
    

    double StarX = 0;
    double StarY = 0;
    double CurrX = 0;
    double CurrY = 0; //delicious
    double NextX = 0;
    double NextY = 0;
    double RealX = 0; //Necessary to revert back to the original coordinate system
    const int nSteps = 25000;
    double step = 0.01;
    double time = 2*60*60;
    
    double stepAzList[nSteps];
    double stepZdList[nSteps];
        
    for(unsigned long int i=0; i<predictionStartPointsAzList.size(); i++){
	   	StarX = predictionStartPointsZdList[i] * TMath::Cos(predictionStartPointsAzList[i]) - CelestialPoleX;
    	StarY = predictionStartPointsZdList[i] * TMath::Sin(predictionStartPointsAzList[i]) - CelestialPoleY;
		CurrX = StarX;
		CurrY = StarY;
		for(int j = 0; j < nSteps; ++j){
			if(CurrX > 0 && CurrY > 0){	//Quad 1
				NextX = CurrX - step;
				NextY = TMath::Sqrt(std::pow(CurrX,2) + std::pow(CurrY,2) - std::pow(NextX,2));
				RealX = NextX + CelestialPoleX; //Remember to change if ever the coordinate system is modified
                CurrX = NextX;
                CurrY = NextY;
                stepZdList[j] = TMath::Sqrt(RealX*RealX + CurrY*CurrY);
                stepAzList[j] = std::atan2(CurrY, RealX) + directionOfNorth * degToRad;	
			}
			if(CurrX < 0 && CurrY > 0){	//Quad 2
				NextY = CurrY - step;
				NextX = -TMath::Sqrt(std::pow(CurrX,2) + std::pow(CurrY,2) - std::pow(NextY,2));
				RealX = NextX + CelestialPoleX; //Remember to change if ever the coordinate system is modified
                CurrX = NextX;
                CurrY = NextY;
                stepZdList[j] = TMath::Sqrt(RealX*RealX + CurrY*CurrY);
                stepAzList[j] = std::atan2(CurrY, RealX) + directionOfNorth * degToRad;	
			}
			if(CurrX < 0 && CurrY < 0){	//Quad 3
				NextX = CurrX + step;
				NextY = -TMath::Sqrt(std::pow(CurrX,2) + std::pow(CurrY,2) - std::pow(NextX,2));
				RealX = NextX + CelestialPoleX; //Remember to change if ever the coordinate system is modified
                CurrX = NextX;
                CurrY = NextY;
                stepZdList[j] = TMath::Sqrt(RealX*RealX + CurrY*CurrY);
                stepAzList[j] = std::atan2(CurrY, RealX) + directionOfNorth * degToRad;
			}
			if(CurrX > 0 && CurrY < 0){	//Quad 4
				NextY = CurrY + step;
				NextX = TMath::Sqrt(std::pow(CurrX,2) + std::pow(CurrY,2) - std::pow(NextY,2));
				RealX = NextX + CelestialPoleX; //Remember to change if ever the coordinate system is modified
                CurrX = NextX;
                CurrY = NextY;
                stepZdList[j] = TMath::Sqrt(RealX*RealX + CurrY*CurrY);
                stepAzList[j] = std::atan2(CurrY, RealX) + directionOfNorth * degToRad;
			}
		}
		TGraphPolar *circle = new TGraphPolar(nSteps, stepAzList, stepZdList);
        circle->SetPolargram(polarFrame);
        circle->SetLineColor(kBlue);
        circle->Draw("L SAME"); // "L" = line connecting points
        
        
        
    }
    
    
    
    
    
    
    
    //ROOT graphing tool. THESE LINES MUST BE AT THE END OF THE PAGE ALWAYS.
    
    c1.Update();
    app.Run();

    return 0;
}
