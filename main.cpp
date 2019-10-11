#include "mkl.h"
#include "pzmatrix.h"
#include "TPZPardisoControl.h"
#include "pzcmesh.h"
#include "pzgengrid.h"
#include "TPZMatLaplacian.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include "TPZVTKGeoMesh.h"
#include "pzgnode.h"
#include "TPZRefPattern.h"
#include "TPZRefPatternTools.h"
#include "TPZGeoElRefPattern.h"
#include "pzgmesh.h"

enum MMATID {EGroup, Ebc1, Ebc2, Ebc3, Ebc4, Emat1};

int main(){
    
    MElementType element = EQuadrilateral;
    bool solve = true;
    bool RefinementPattern = false;
    int porder = 1;
    
    if(element == EQuadrilateral){
        
        TPZGeoMesh *gmesh = new TPZGeoMesh();
        gmesh->NodeVec().Resize(6);
        
        for (int i=0; i<6; i++) {
            TPZManVector<REAL,4> x(3,0);
            if (i>2) {
                x[1]=1;
            }
            if (i==1 || i==4){
                x[0]=1;
            } else if (i==2 || i==5){
                x[0]=2;
            }
            
            TPZGeoNode node;
            node.Initialize(x, *gmesh);
            gmesh->NodeVec()[i] = node;
//            gmesh->NodeVec()[i].SetCoord(i,x);
//            gmesh->NodeVec()[i].SetNodeId(i);
        }
        
        TPZVec<int64_t> cornerindexes(4,0);
        cornerindexes[0]=0;cornerindexes[1]=1;cornerindexes[2]=4;cornerindexes[3]=3;
        int64_t index = 0;
        TPZGeoEl *gel = gmesh->CreateGeoElement(EQuadrilateral, cornerindexes, EGroup, index);
        
        cornerindexes[0]=1;cornerindexes[1]=2;cornerindexes[2]=5;cornerindexes[3]=4;
        index=1;
        gel = gmesh->CreateGeoElement(EQuadrilateral, cornerindexes, EGroup, index);
        
        gmesh->BuildConnectivity();
        TPZVec<TPZGeoEl *> subel;
        gmesh->Element(0)->Divide(subel);
        gmesh->BuildConnectivity();
        
        gmesh->Print();
        std::ofstream out("Geometry.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(gmesh, out, true);
        
        if(RefinementPattern){
            TPZGeoMesh *gmesh2 = new TPZGeoMesh(*gmesh);
            gmesh2->BuildConnectivity();
            
            int nnodes = gmesh->NNodes();
            gmesh->NodeVec().Resize(nnodes +2);
            TPZVec<REAL> x(3,0);
            x[0]=-0.5;x[1]=-1;
            gmesh->NodeVec()[nnodes].Initialize(nnodes, x, *gmesh);
            x[0]=0.9;x[1]=1;
            gmesh->NodeVec()[nnodes+1].Initialize(nnodes+1, x, *gmesh);
            
            cornerindexes[0]=0;cornerindexes[1]=1;cornerindexes[2]=4;cornerindexes[3]=5;
            index = 1;
            
            gmesh->CreateGeoElement(EQuadrilateral, cornerindexes, EGroup, index);
            cornerindexes[0]=1;cornerindexes[1]=2;cornerindexes[2]=3;cornerindexes[3]=4;
            index = 2;
            gmesh->CreateGeoElement(EQuadrilateral, cornerindexes, EGroup, index);
            
            gmesh->Element(1)->SetFather(gmesh->Element(0));
            gmesh->Element(2)->SetFather(gmesh->Element(0));
            
            gmesh->BuildConnectivity();
            gmesh->Print(std::cout);
            std::ofstream out("Geometry.vtk");
            TPZVTKGeoMesh vtk;
            vtk.PrintGMeshVTK(gmesh, out, true);
            
            TPZRefPattern *ref = new TPZRefPattern(*gmesh);
            ref->SetName("Karol");
            TPZAutoPointer<TPZRefPattern> autoref(ref);
            gRefDBase.InsertRefPattern(autoref);
            autoref->InsertPermuted();
            auto karol = gRefDBase.FindRefPattern("Karol");
            
            TPZVec<TPZGeoEl *> gelvec;
//            gmesh2->Element(0)->SetRefPattern(karol);
            gmesh2->Element(0)->Divide(gelvec);
            gmesh2->Element(1)->SetRefPattern(karol);
            gmesh2->Element(1)->Divide(gelvec);
            std::ofstream out2("Geometry2.vtk");
            TPZVTKGeoMesh vtk2;
            vtk2.PrintGMeshVTK(gmesh2, out2, true);
            
            for(int i=4;i<6;i++){
                std::map<int, std::pair<TPZGeoEl *, std::map<int,int> > > neighCorresp;
                TPZRefPatternTools::ModelRefPattern(gmesh2->Element(i),neighCorresp);
                TPZAutoPointer<TPZRefPattern> ref = TPZRefPatternTools::DragModelPatNodes(gmesh2->Element(i), karol,neighCorresp);
                gmesh2->Element(i)->SetRefPattern(ref);
//                gmesh2->Element(1)->SetRefPattern();
                gmesh2->Element(i)->Divide(gelvec);
                std::ofstream out2("Geometry2.vtk");
                TPZVTKGeoMesh vtk2;
                vtk2.PrintGMeshVTK(gmesh2, out2, true);
            }
            
//            for(int i=0;i<600;i++){
//                gmesh2->Element(i)->SetRefPattern(karol);
//                gmesh2->Element(i)->Divide(gelvec);
//            }
            std::cout << gmesh2->NElements() << std::endl;
            std::ofstream out3("Geometry3.vtk");
            TPZVTKGeoMesh vtk3;
            vtk3.PrintGMeshVTK(gmesh2, out3, true);
        }
        
        if(solve){
            gmesh->Element(0)->CreateBCGeoEl(4, -Ebc1);
            gmesh->Element(0)->CreateBCGeoEl(6, -Ebc1);
            gmesh->Element(0)->CreateBCGeoEl(7, -Ebc1);
            gmesh->Element(1)->CreateBCGeoEl(4, -Ebc1);
            gmesh->Element(1)->CreateBCGeoEl(5, -Ebc1);
            gmesh->Element(1)->CreateBCGeoEl(6, -Ebc1);

            TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
            cmesh->SetDefaultOrder(porder);
            cmesh->SetDimModel(2);
            
            TPZMatLaplacian *matloc = new TPZMatLaplacian(EGroup);
            matloc->SetDimension(2);
            matloc->SetSymmetric();
            
            TPZMaterial *material;
            material = matloc;
            int nstate = 1;
            
            TPZFMatrix<STATE> val1(nstate,nstate,0.), val2(nstate,1,0.);
            
            TPZMaterial * BCond1 = material->CreateBC(material,-Ebc1,0, val1, val2);
            
            cmesh->InsertMaterialObject(material);
            cmesh->InsertMaterialObject(BCond1);
//            cmesh->SetAllCreateFunctionsContinuous();
            cmesh->SetAllCreateFunctionsHDiv();
            cmesh->AutoBuild();
            
            std::ofstream file("cmeshhdiv.dat");
            cmesh->Print(file);
            
            TPZAnalysis * Analysis = new TPZAnalysis(cmesh,true);
            Analysis->Run();
            
            TPZFMatrix<REAL> sol = Analysis->Solution();
            
            bool plotshape = true;
            if(plotshape)
            {
                TPZFMatrix<REAL> sol0 = sol;
                for (int i=0; i<sol0.Rows() ;i++){
                    
                    TPZFNMatrix<3,REAL> sol = cmesh->Solution();
                    sol.Zero();
                    sol(i,0) = 1;
                    
                    cmesh->LoadSolution(sol);
                    Analysis->LoadSolution(sol);
                    
                    TPZStack<std::string> vecnames,scalnames;
                    // scalar
                    scalnames.Push("State");
                    Analysis->DefineGraphMesh(2, scalnames, vecnames, "../ShapeFunctions.vtk");
                    Analysis->PostProcess(3);
                }
                cmesh->LoadSolution(sol0);
                Analysis->LoadSolution(sol0);
            }
        }
        
    } else if (element == ETriangle){
        TPZManVector<REAL,4> x0(3,-1.),x1(3,1.),x2(3,1.);
        x0[0] = 0;
        x0[1] = 0;
        x0[2] = 0;
        
        x1[0] = 1;
        x1[1] = 0;
        x1[2] = 0;
        
        x2[0] = 0;
        x2[1] = 1;
        x2[2] = 0;
        
        TPZGeoMesh *gmesh = new TPZGeoMesh();
        gmesh->NodeVec().Resize(3);
        
        gmesh->NodeVec()[0].SetCoord(x0);
        gmesh->NodeVec()[1].SetCoord(x1);
        gmesh->NodeVec()[2].SetCoord(x2);
        
        gmesh->NodeVec()[0].SetNodeId(0);
        gmesh->NodeVec()[1].SetNodeId(1);
        gmesh->NodeVec()[2].SetNodeId(2);
        
        TPZVec<int64_t> cornerindexes(3,0);
        cornerindexes[0]=0;cornerindexes[1]=1;cornerindexes[2]=2;
        
        int64_t index = 0;
        gmesh->CreateGeoElement(ETriangle, cornerindexes, EGroup, index);
        
        gmesh->BuildConnectivity();
        gmesh->Element(0)->CreateBCGeoEl(3, -Ebc1);
        gmesh->Element(0)->CreateBCGeoEl(4, -Ebc2);
        gmesh->Element(0)->CreateBCGeoEl(5, -Ebc3);
        
        std::ofstream out("Geometry.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(gmesh, out, true);
        
        TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
        cmesh->SetDefaultOrder(porder);
        
        TPZMatLaplacian *matloc = new TPZMatLaplacian(EGroup);
        matloc->SetDimension(2);
        matloc->SetSymmetric();
        
        TPZMaterial *material;
        material = matloc;
        int nstate = 1;
        
        TPZFMatrix<STATE> val1(nstate,nstate,0.), val2(nstate,1,0.);
        
        TPZMaterial * BCond1 = material->CreateBC(material,-Ebc1,0, val1, val2);
        TPZMaterial * BCond2 = material->CreateBC(material,-Ebc2,0, val1, val2);
        TPZMaterial * BCond3 = material->CreateBC(material,-Ebc3,0, val1, val2);
        
        cmesh->InsertMaterialObject(material);
        cmesh->InsertMaterialObject(BCond1);
        cmesh->InsertMaterialObject(BCond2);
        cmesh->InsertMaterialObject(BCond3);
        cmesh->SetAllCreateFunctionsContinuous();
        cmesh->AutoBuild();
        
        TPZAnalysis * Analysis = new TPZAnalysis(cmesh,true);
        Analysis->Run();
        
        TPZFMatrix<REAL> sol = Analysis->Solution();
        
        bool plotshape = true;
        if(plotshape)
        {
            TPZFMatrix<REAL> sol0 = sol;
            for (int i=0; i<sol0.Rows() ;i++){
                
                TPZFNMatrix<3,REAL> sol = cmesh->Solution();
                sol.Zero();
                sol(i,0) = 1;
                
                cmesh->LoadSolution(sol);
                Analysis->LoadSolution(sol);
                
                TPZStack<std::string> vecnames,scalnames;
                // scalar
                scalnames.Push("State");
                Analysis->DefineGraphMesh(2, scalnames, vecnames, "../ShapeFunctions.vtk");
                Analysis->PostProcess(3);
            }
            cmesh->LoadSolution(sol0);
            Analysis->LoadSolution(sol0);
        }
        
    }
    
    return 0;
}
