#include <iostream>
#include <fstream>
#include <TGraphPolar.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TGraphPolargram.h>
#include <TMath.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TEllipse.h>


std::string LabelsPolarTextual[24] = {
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

    auto pi = TMath::Pi();
    double toRad = pi / 180;

    TApplication app("app", 0,0);
    
    gStyle->SetOptTitle(0);
    TCanvas c1;
    auto polarFrame = new TGraphPolargram("polarFrame", 0, 360*toRad, 0, 90);
    polarFrame->SetToRadian();
    polarFrame->SetNdivRadial(9);
    polarFrame->SetNdivPolar(24);
    polarFrame->Draw("P SAME");
    polarFrame->SetPolarLabelSize(0);

    double r = 1.1; // Relative radius (between 0 and 1)
    double base_angle = 90;

    int PolarDivs = polarFrame->GetNdivPolar();
    int RadialDivs = polarFrame->GetNdivRadial();
    double angle_deg = 360.0/PolarDivs; // Each div

    for(int i = 0; i<PolarDivs; i++){
        double x = r * TMath::Cos((base_angle + i * angle_deg) * TMath::Pi() / 180.0);
        double y = r * TMath::Sin((base_angle + i * angle_deg) * TMath::Pi() / 180.0);
        TLatex *Label = new TLatex(x, y, LabelsPolarTextual[i].c_str());
        Label->SetTextAlign(22); // Center align
        Label->SetTextSize(0.05);
        Label->Draw();
    }

    for(int i = 0; i<RadialDivs; i++){
        double y = (i+1) * 0.1 + 0.0125*i;
        double x = 0;
        TLatex *Label = new TLatex(x, y, std::to_string((i+1)*10).c_str());
        Label->SetTextAlign(25); // Center align
        Label->SetTextSize(0.02);
        Label->Draw();
    }

    std::ifstream file("taken.dat");
    std::ifstream file2("predicted.dat");


    double az = 0, zd = 0;
    std::vector<double> measaz;
    std::vector<double> measzd;
    
    while(file >> az >> zd){
        measaz.push_back((az + 90) * toRad);
        measzd.push_back(zd);
    }
    
    //TAKEN STARS
    Double_t azGraph[measaz.size()];
    Double_t zdGraph[measzd.size()];

    for(unsigned long int i=0; i<measaz.size(); i++){
        azGraph[i] = measaz[i];
        zdGraph[i] = measzd[i];
    }
    TGraphPolar gr1(measaz.size(), azGraph, zdGraph);

    gr1.SetPolargram(polarFrame);
    gr1.SetMarkerStyle(29);
    gr1.SetMarkerColor(kBlack);
    gr1.SetMinRadial(0);
    gr1.SetMaxRadial(90);
    gr1.SetMarkerSize(1.5);
    gr1.Draw("P SAME");

    //PREDICTION
    
    double azInit = 0, zdInit = 0, azFinal = 0, zdFinal = 0;
    
    std::vector<double> StartPointsAz;
    std::vector<double> EndPointsAz;
    std::vector<double> StartPointsZd;
    std::vector<double> EndPointsZd;
    std::vector<double> AzProjAll;
    std::vector<double> ZdProjAll;
    
    
    while(file2 >> azInit >> zdInit >> azFinal >> zdFinal ){
        StartPointsAz.push_back((azInit+90) *toRad);
        StartPointsZd.push_back(90 - zdInit);
        EndPointsAz.push_back((azFinal+90) *toRad);
        EndPointsZd.push_back(90 - zdFinal);
        AzProjAll.push_back((azInit+90) *toRad);
        AzProjAll.push_back((azFinal+90) *toRad);
        ZdProjAll.push_back(90 - zdInit);
        ZdProjAll.push_back(90 - zdFinal);
    }

    Double_t AzProjGraph[AzProjAll.size()];
    Double_t ZdProjGraph[ZdProjAll.size()];

    for(unsigned long int i=0; i<AzProjAll.size(); i++){
        AzProjGraph[i] = AzProjAll[i];
        ZdProjGraph[i] = ZdProjAll[i];
    }


    TGraphPolar gr2(AzProjAll.size(), AzProjGraph, ZdProjGraph);

    gr2.SetPolargram(polarFrame);
    gr2.SetMarkerStyle(29);
    gr2.SetMarkerColor(kRed);
    gr2.SetMinRadial(0);
    gr2.SetMaxRadial(90);
    gr2.SetMarkerSize(1.5);
    gr2.Draw("P SAME");

    double CelestialPoleX = 0;
    double CelestialPoleY = 44.5892595; //Latitude of the observatory :)
    double StarX = 0;
    double StarY = 0;

    for (size_t i = 0; i < StartPointsAz.size(); ++i) {
        StarX = -( StartPointsZd[i] * TMath::Cos(StartPointsAz[i]) - CelestialPoleX );
        StarY = StartPointsZd[i] * TMath::Sin(StartPointsAz[i]) - CelestialPoleY;

        double CurrX = StarX;
        double CurrY = StarY;
        double NextX = CurrX;
        double NextY = CurrY;
        double RealY;
        const int nPoints = 25000;
        double step = 0.001;
        
        double angle[nPoints];
        double radius[nPoints];
        for(int j = 0; j < nPoints; ++j){
            if(CurrX > 0 && CurrY > 0){
                NextY = CurrY - step;
                NextX = TMath::Sqrt(CurrX*CurrX + CurrY*CurrY - NextY*NextY);
                CurrX = NextX;
                CurrY = NextY;
                RealY = CurrY + CelestialPoleY;

                radius[j] = TMath::Sqrt(CurrX*CurrX + RealY*RealY);
                angle[j] = std::atan2(CurrX, RealY) + 90*toRad;
            }
            if(CurrX < 0 && CurrY > 0){
                NextX = CurrX + step;
                NextY = TMath::Sqrt(CurrY*CurrY + CurrX*CurrX - NextX*NextX);
                CurrX = NextX;
                CurrY = NextY;
                RealY = CurrY + CelestialPoleY;

                radius[j] = TMath::Sqrt(CurrX*CurrX + RealY*RealY);
                angle[j] = std::atan2(CurrX, RealY) + 90*toRad;
            }
            if(CurrX > 0 && CurrY < 0){
                NextX = CurrX - step;
                NextY = -TMath::Sqrt(CurrY*CurrY + CurrX*CurrX - NextX*NextX);
                CurrX = NextX;
                CurrY = NextY;
                RealY = CurrY + CelestialPoleY;

                radius[j] = TMath::Sqrt(CurrX*CurrX + RealY*RealY);
                angle[j] = std::atan2(CurrX, RealY) + 90*toRad;
            }
            if(CurrX < 0 && CurrY < 0){
                NextY = CurrY + step;
                NextX = - TMath::Sqrt(CurrY*CurrY + CurrX*CurrX - NextY*NextY);
                CurrX = NextX;
                CurrY = NextY;
                RealY = CurrY + CelestialPoleY;

                radius[j] = TMath::Sqrt(CurrX*CurrX + RealY*RealY);
                angle[j] = std::atan2(CurrX, RealY) + 90*toRad;
            }
        }
        TGraphPolar *circle = new TGraphPolar(nPoints, angle, radius);
        circle->SetPolargram(polarFrame);
        circle->SetLineColor(kRed);
        circle->Draw("L SAME"); // "L" = line connecting points

    }

    
    c1.Update();
    app.Run();

    return 0;
}
