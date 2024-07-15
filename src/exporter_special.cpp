#include "exporter.h"
//#include "SubdivEvaluator.h"
//#include "adolcsubdiv.h"

#include "Common.h"
#include <fbxsdk.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>

#define SCALAR (1.0)
#define VE

static FbxNode *CreateCubeMesh(FbxScene *pScene, const char *pName);
static FbxNode *createBaseMesh(FbxScene *pScene, const char *pName,
                               Model *model, int personId);
static FbxNode *createJointMesh(FbxScene *pScene, const char *pName,
                                FbxAnimLayer *myAnimBaseLayer, Model *model,
                                int personId);
static FbxNode *createSkeleton(FbxScene *pScene, const char *pName,
                               FbxAnimLayer *myAnimBaseLayer, int *mode,
                               FbxNode *pMeshNode, Model *model,
                               FbxNode **bones, int personId, int poseId,
                               int endId);
static void LinkPatchToSkeleton(FbxScene *pScene, FbxNode *pPatch,
                                FbxNode *pSkeletonRoot, Model *model,
                                FbxNode **bones);

static bool convertFbx(std::string path, std::string outpath) {
  FbxManager *lSdkManager = NULL;
  FbxScene *lScene = NULL;
  InitializeSdkObjects(lSdkManager, lScene);

  LoadScene(lSdkManager, lScene, path.c_str());
  SaveScene(lSdkManager, lScene, outpath.c_str(), -1);

  DestroySdkObjects(lSdkManager, true);
  return true;
}

bool writePersonFbx(std::string path, int *mode, Model *model, int personId,
                    int poseId, int endId) {

  FbxManager *lSdkManager = NULL;
  FbxScene *lScene = NULL;
  FbxString pCubeName = "test";
  FbxNode *bones[model->skeleton->getJointCount()];
  // Prepare the FBX SDK.
  InitializeSdkObjects(lSdkManager, lScene);
  // lScene->GetGlobalSettings().SetSystemUnit(FbxSystemUnit::m);

  FbxNode *base = FbxNode::Create(lScene, "base");
  // base->LclScaling.Set(FbxVector4(100, 100, 100));
  lScene->GetRootNode()->AddChild(base);

  FbxNode *lCube = createBaseMesh(lScene, "test", model, personId);
  base->AddChild(lCube);

  FbxAnimStack *myAnimStack = FbxAnimStack::Create(lScene, "My stack");
  FbxAnimLayer *myAnimBaseLayer = FbxAnimLayer::Create(lScene, "Layer0");
  myAnimStack->AddMember(myAnimBaseLayer);
  FbxNode *lJoint =
      createJointMesh(lScene, "joints", myAnimBaseLayer, model, personId);
  base->AddChild(lJoint);
  FbxNode *skel = createSkeleton(lScene, "", myAnimBaseLayer, mode, lCube,
                                 model, bones, personId, poseId, endId);
  base->AddChild(skel);

  LinkPatchToSkeleton(lScene, lCube, bones[0], model, bones);
  // base->LclScaling.Set(FbxVector4(100, 100, 100));

  std::string outpath = path + ".fbx";
  // std::string outpath = path + (char)('X' + mode[0]) + (char)('X' + mode[1])
  // + (char)('X' + mode[2]) + ".fbx"; std::string outpath2 = path + "bin" +
  // (char)('X' + mode[0]) + (char)('X' + mode[1]) + (char)('X' + mode[2]) +
  // ".fbx"; SaveScene(lSdkManager, lScene, outpath.c_str(), -1);
  SaveScene(lSdkManager, lScene, outpath.c_str(), 0);

  DestroySdkObjects(lSdkManager, true);

  return true;
}

static void LinkPatchToSkeleton(FbxScene *pScene, FbxNode *pPatch,
                                FbxNode *pSkeletonRoot, Model *model,
                                FbxNode **bones) {
  auto skel = model->skeleton;

  int i, j;
  FbxAMatrix lXMatrix;
  FbxGeometry *lPatchAttribute = (FbxGeometry *)pPatch->GetNodeAttribute();
  FbxSkin *lSkin = FbxSkin::Create(pScene, "");
  lPatchAttribute->AddDeformer(lSkin);

  for (int i = 0; i < skel->getJointCount(); i++) {

    FbxCluster *lClusterToRoot = FbxCluster::Create(pScene, "");
    lClusterToRoot->SetLink(bones[i]);
    lClusterToRoot->SetLinkMode(FbxCluster::eTotalOne);
    for (int j = 0; j < model->mean.cols(); j++) {
      for (int k = 0; k < 4; k++) {
        int influenceJoint = model->closestJoints(k, j);
        //#ifdef POSE_DEFORM
        double weight = model->weights(influenceJoint, j);
        // weight = .25;
        //#else
        //				double weight = global->initialSkinWeights(k,
        //j); #endif
        if (influenceJoint != i)
          continue;
        lClusterToRoot->AddControlPointIndex(j, weight);
      }

      /*int influenceJoint = global->closestJoints(0, j);
      if (influenceJoint != i) continue;
      lClusterToRoot->AddControlPointIndex(j, 1);*/
    }

    FbxScene *lScene = pPatch->GetScene();
    // lXMatrix = pPatch->EvaluateGlobalTransform();
    // lClusterToRoot->SetTransformMatrix(lXMatrix);
    // std::cout << lXMatrix << std::endl;
    lXMatrix = bones[i]->EvaluateGlobalTransform();
    lClusterToRoot->SetTransformLinkMatrix(lXMatrix);
    lSkin->AddCluster(lClusterToRoot);
  }
}

static FbxNode *createSkeleton(FbxScene *pScene, const char *pName,
                               FbxAnimLayer *myAnimBaseLayer, int *mode,
                               FbxNode *pMeshNode, Model *model,
                               FbxNode **bones, int personId, int poseId,
                               int endId) {
  auto skel = model->skeleton;
  Matrix3X jointPos(3, skel->getJointCount());
  jointPos.setZero();
  Person &person = model->persons[personId];

  jointPos = model->joints;

  for (int j = 0; j < model->pcs.size(); j++) {
    jointPos += person.coeff[j] * model->jointpcs[j];
  }

  // std::cout << "max: " << jointPos.maxCoeff() << " min: " <<
  // jointPos.minCoeff() << std::endl;

  FbxSkeleton *bone = FbxSkeleton::Create(pScene, pName);
  bone->SetSkeletonType(FbxSkeleton::eRoot);
  FbxNode *start = FbxNode::Create(pScene, "proxy");
  start->SetNodeAttribute(bone);

  for (int i = 0; i < skel->getJointCount(); i++) {
    auto pair = skel->getJoint(i);
    FbxString lRootName("Skeleton_");
    FbxString partName(pair.second.c_str());
    lRootName += partName;
    FbxSkeleton *bone = FbxSkeleton::Create(pScene, pName);
    // if (i == 0) {
    //	bone->SetSkeletonType(FbxSkeleton::eRoot);
    //}
    /*
    else if (i == 1 || i == 4 || i == 7 || i == 12 || i == 15) {
    bone->SetSkeletonType(FbxSkeleton::eEffector);
    }
    */
    bone->SetSkeletonType(FbxSkeleton::eLimbNode);

    bone->LimbLength.Set(0.001);
    bones[i] = FbxNode::Create(pScene, lRootName.Buffer());
    bones[i]->SetNodeAttribute(bone);
    bones[i]->SetRotationOrder(FbxNode::EPivotSet::eDestinationPivot,
                               FbxEuler::eOrderXYZ);
  }

  bones[0]->LclTranslation.Set(
      FbxVector4(jointPos(0, 0), jointPos(1, 0), jointPos(2, 0)));

  for (int i = 1; i < skel->getJointCount(); i++) {
    auto pair = skel->getJoint(i);
    int parent = skel->getParentJoint(i);
    bones[parent]->AddChild(bones[i]);
    Vector3 trans = jointPos.col(i) - jointPos.col(parent);
    // std::cout << i << ":" << parent << std::endl;
    bones[i]->LclTranslation.Set(FbxVector4(trans[0], trans[1], trans[2]));
  }

  start->AddChild(bones[0]);

  {

    // auto pair = skel->getJoint(i);
    FbxAnimCurveNode *myAnimCurveNode =
        bones[0]->LclTranslation.GetCurveNode(myAnimBaseLayer, true);

    FbxAnimCurve *myTransCurveX = bones[0]->LclTranslation.GetCurve(
        myAnimBaseLayer, FBXSDK_CURVENODE_COMPONENT_X, true);
    FbxAnimCurve *myTransCurveY = bones[0]->LclTranslation.GetCurve(
        myAnimBaseLayer, FBXSDK_CURVENODE_COMPONENT_Y, true);
    FbxAnimCurve *myTransCurveZ = bones[0]->LclTranslation.GetCurve(
        myAnimBaseLayer, FBXSDK_CURVENODE_COMPONENT_Z, true);

    myTransCurveX->KeyModifyBegin();
    myTransCurveY->KeyModifyBegin();
    myTransCurveZ->KeyModifyBegin();

    FbxTime myTime;
    int myKeyIndex = -1;
    myTime.SetSecondDouble((myKeyIndex + 1) * 0.040);
    myKeyIndex = myTransCurveX->KeyAdd(myTime);
    myTransCurveY->KeyAdd(myTime);
    myTransCurveZ->KeyAdd(myTime);
    myTransCurveX->KeySetValue(myKeyIndex, 0);
    myTransCurveY->KeySetValue(myKeyIndex, 0);
    myTransCurveZ->KeySetValue(myKeyIndex, 0);
    myTransCurveX->KeySetInterpolation(myKeyIndex,
                                       FbxAnimCurveDef::eInterpolationLinear);
    myTransCurveY->KeySetInterpolation(myKeyIndex,
                                       FbxAnimCurveDef::eInterpolationLinear);
    myTransCurveZ->KeySetInterpolation(myKeyIndex,
                                       FbxAnimCurveDef::eInterpolationLinear);

    for (int t = poseId; t < endId;
         t++) { // fix person index model->persons[s]->scans.size()
      MatrixX trans = model->scans[t].trans;

      myTime.SetSecondDouble((myKeyIndex + 1) * 0.04);
      myKeyIndex = myTransCurveX->KeyAdd(myTime);
      myTransCurveY->KeyAdd(myTime);
      myTransCurveZ->KeyAdd(myTime);

      myTransCurveX->KeySetValue(myKeyIndex, trans(0, 0));
      myTransCurveY->KeySetValue(myKeyIndex, trans(1, 0));
      myTransCurveZ->KeySetValue(myKeyIndex, trans(2, 0));

      myTransCurveX->KeySetInterpolation(myKeyIndex,
                                         FbxAnimCurveDef::eInterpolationLinear);
      myTransCurveY->KeySetInterpolation(myKeyIndex,
                                         FbxAnimCurveDef::eInterpolationLinear);
      myTransCurveZ->KeySetInterpolation(myKeyIndex,
                                         FbxAnimCurveDef::eInterpolationLinear);
    }

    myTransCurveX->KeyModifyEnd();
    myTransCurveY->KeyModifyEnd();
    myTransCurveZ->KeyModifyEnd();
  }

  for (int i = 0; i < skel->getJointCount(); i++) {
    auto pair = skel->getJoint(i);
    FbxAnimCurveNode *myAnimCurveNode =
        bones[i]->LclRotation.GetCurveNode(myAnimBaseLayer, true);

    FbxAnimCurve *myRotCurveX = bones[i]->LclRotation.GetCurve(
        myAnimBaseLayer, FBXSDK_CURVENODE_COMPONENT_X, true);
    FbxAnimCurve *myRotCurveY = bones[i]->LclRotation.GetCurve(
        myAnimBaseLayer, FBXSDK_CURVENODE_COMPONENT_Y, true);
    FbxAnimCurve *myRotCurveZ = bones[i]->LclRotation.GetCurve(
        myAnimBaseLayer, FBXSDK_CURVENODE_COMPONENT_Z, true);

    myRotCurveX->KeyModifyBegin();
    myRotCurveY->KeyModifyBegin();
    myRotCurveZ->KeyModifyBegin();

    FbxTime myTime;
    myTime.SetSecondDouble(0.0);
    int myKeyIndex = -1;

    myTime.SetSecondDouble((myKeyIndex + 1) * 0.04);
    myKeyIndex = myRotCurveX->KeyAdd(myTime);
    myRotCurveY->KeyAdd(myTime);
    myRotCurveZ->KeyAdd(myTime);

    myRotCurveX->KeySetValue(myKeyIndex, 0);
    myRotCurveY->KeySetValue(myKeyIndex, 0);
    myRotCurveZ->KeySetValue(myKeyIndex, 0);

    myRotCurveX->KeySetInterpolation(myKeyIndex,
                                     FbxAnimCurveDef::eInterpolationLinear);
    myRotCurveY->KeySetInterpolation(myKeyIndex,
                                     FbxAnimCurveDef::eInterpolationLinear);
    myRotCurveZ->KeySetInterpolation(myKeyIndex,
                                     FbxAnimCurveDef::eInterpolationLinear);

    for (int t = poseId; t < endId;
         t++) { // fix scan size global->persons[s]->scans.size()
      /*
                  Matrix3 rot = model->scans[t].posebase.col(i);
                  if (t == 0) {
                          std::cout << i << " : " <<
         model->skeleton->getJoint(i).second << " -> " <<
         model->skeleton->getJoint(i).first << " : \n"; std::cout <<  ":" << rot
         << std::endl;
                  }
                  //rot.transposeInPlace();
                  auto euler = 180.0 / 3.141 * rot.eulerAngles(mode[0], mode[1],
         mode[2]);
      */
      auto euler = 180.0 / 3.141 * model->scans[t].posebase.col(i);
      myTime.SetSecondDouble((myKeyIndex + 1) * 0.04);
      myKeyIndex = myRotCurveX->KeyAdd(myTime);
      myRotCurveY->KeyAdd(myTime);
      myRotCurveZ->KeyAdd(myTime);

      myRotCurveX->KeySetValue(myKeyIndex, euler[2]);
      myRotCurveY->KeySetValue(myKeyIndex, euler[1]);
      myRotCurveZ->KeySetValue(myKeyIndex, euler[0]);
      /*
      std::cout << s << " : " << t << " : " << pair.first << " -> " <<
              pair.second << " : " <<
              std::fixed << std::setprecision(2) <<
              euler[0] << ", " <<
              std::fixed << std::setprecision(2) <<
              euler[1] << ", " <<
              std::fixed << std::setprecision(2) <<
              euler[2] << std::endl;
      */
      myRotCurveX->KeySetInterpolation(myKeyIndex,
                                       FbxAnimCurveDef::eInterpolationLinear);
      myRotCurveY->KeySetInterpolation(myKeyIndex,
                                       FbxAnimCurveDef::eInterpolationLinear);
      myRotCurveZ->KeySetInterpolation(myKeyIndex,
                                       FbxAnimCurveDef::eInterpolationLinear);
    }

    myRotCurveX->KeyModifyEnd();
    myRotCurveY->KeyModifyEnd();
    myRotCurveZ->KeyModifyEnd();
  }

  return start;
}

FbxNode *createJointMesh(FbxScene *pScene, const char *pName,
                         FbxAnimLayer *myAnimBaseLayer, Model *model,
                         int personId) {
  FbxMesh *lMesh = FbxMesh::Create(pScene, pName);

  FbxNode *lNode = FbxNode::Create(pScene, pName);

  // std::cout << "max: " << poseRegressor.maxCoeff() << " min: " <<
  // poseRegressor.minCoeff() << std::endl;

  Matrix3X a;
  a = model->joints;
  Person &person = model->persons[personId];

  // for (int j = 0; j < model->pcs.size(); j++) {
  //	a += person.coeff[j] * model->pcs[j];
  //}

  // Eigen::Matrix<int, 5, Eigen::Dynamic> f = model->top.quads;
  lMesh->InitControlPoints(a.cols());
  FbxVector4 *lControlPoints = lMesh->GetControlPoints();

  for (int i = 0; i < model->joints.cols(); i++) {
    lControlPoints[i] = FbxVector4(a(0, i), a(1, i), a(2, i));
  }
  /*
  for (int i = 0; i < f.cols(); i++)
  {
          lMesh->BeginPolygon(-1, -1, -1, false);

          for (int j = 0; j < 4; j++)
          {
                  lMesh->AddPolygon(f(j, i));
          }

          lMesh->EndPolygon();
  }*/

  for (int i = 0; i < model->jointpcs.size(); i++) {
    for (int dir = -1; dir <= 1; dir += 2) {
      std::stringstream buffer;
      buffer << "shape" << std::setfill('0') << std::setw(5) << i;
      std::string blendName = buffer.str();
      FbxShape *shape = FbxShape::Create(pScene, blendName.c_str());
      shape->InitControlPoints(a.cols());
      FbxVector4 *pts = shape->GetControlPoints();

      for (int j = 0; j < a.cols(); j++) {
        Vector3 pt = 10 * model->jointpcs[i].col(j);
        pts[j] = FbxVector4(dir * pt[0] + a(0, j), dir * pt[1] + a(1, j),
                            dir * pt[2] + a(2, j));
      }

      FbxBlendShape *lBlendShape = FbxBlendShape::Create(pScene, "");
      FbxBlendShapeChannel *lBlendShapeChannel =
          FbxBlendShapeChannel::Create(pScene, "");
      lMesh->AddDeformer(lBlendShape);
      lBlendShape->AddBlendShapeChannel(lBlendShapeChannel);
      lBlendShapeChannel->AddTargetShape(shape);
      lBlendShapeChannel->DeformPercent.Set(10 * fabs(person.coeff[i]));
    }
  }

  lNode->SetNodeAttribute(lMesh);
  return lNode;
}

/*
FbxNode* createJointMesh(FbxScene* pScene, char* pName, FbxAnimLayer*
myAnimBaseLayer, Model* model) { FbxMesh* lMesh = FbxMesh::Create(pScene,
pName); FbxNode* lNode = FbxNode::Create(pScene, pName);

        Matrix3X jointPos(3, 16);
        jointPos.setZero();
        auto skel = model->skeleton;

        for (int i = 1; i < 16; i++) {
                jointPos(0, i) = global->jointRegressor(3 * (i - 1) + 0,
0).getValue(); jointPos(1, i) = global->jointRegressor(3 * (i - 1) + 1,
0).getValue(); jointPos(2, i) = global->jointRegressor(3 * (i - 1) + 2,
0).getValue();
        }

        FbxShape* lShape = FbxShape::Create(pScene, "mean shape");
        FbxBlendShape* lBlendShape = FbxBlendShape::Create(pScene, "");
        FbxBlendShapeChannel* lBlendShapeChannel =
FbxBlendShapeChannel::Create(pScene, ""); lMesh->AddDeformer(lBlendShape);
        lBlendShape->AddBlendShapeChannel(lBlendShapeChannel);
        lBlendShapeChannel->AddTargetShape(lShape);
        lBlendShapeChannel->DeformPercent = 0;

        lMesh->InitControlPoints(jointPos.cols());
        lShape->InitControlPoints(jointPos.cols());
        FbxVector4* lControlPoints = lMesh->GetControlPoints();
        FbxVector4* lMeanPoints = lShape->GetControlPoints();

        for (int i = 0; i < jointPos.cols(); i++) {
                lControlPoints[i] = FbxVector4(0, 0, 0);
                lMeanPoints[i] = FbxVector4(jointPos(0, i), jointPos(1, i),
jointPos(2, i));
                //lMeanPoints[i] = FbxVector4(0,0,0);
        }



        FbxAnimCurveNode* myAnimCurveNode =
lBlendShapeChannel->DeformPercent.GetCurveNode(myAnimBaseLayer);

        FbxAnimCurve* myScaleCurve =
lBlendShapeChannel->DeformPercent.GetCurve(myAnimBaseLayer, true);
        myScaleCurve->KeyModifyBegin();
        FbxTime myTime;
        myTime.SetSecondDouble(0.0);
        int myKeyIndex = -1;

        myTime.SetSecondDouble((myKeyIndex + 1)*1.0);
        myKeyIndex = myScaleCurve->KeyAdd(myTime);
        myScaleCurve->KeySetValue(myKeyIndex, 100);
        myScaleCurve->KeySetInterpolation(myKeyIndex,
FbxAnimCurveDef::eInterpolationLinear);


        for (int s = 0; s < global->persons.size(); s++) {
                for (int t = 0; t < global->persons[s]->scans.size(); t++) {

                        double scale = global->persons[s]->scale(0,
0).getValue();
                        //rot.transposeInPlace();

                        myTime.SetSecondDouble((myKeyIndex + 1)*1.0);
                        myKeyIndex = myScaleCurve->KeyAdd(myTime);

                        myScaleCurve->KeySetValue(myKeyIndex, 100 *
SCALAR*(scale));



                        myScaleCurve->KeySetInterpolation(myKeyIndex,
FbxAnimCurveDef::eInterpolationLinear);
                }
        }

        myScaleCurve->KeyModifyEnd();

        for (int i = 0; i < global->component_count; i++) {
                FbxBlendShape* lBlendShape = FbxBlendShape::Create(pScene, "");
                FbxBlendShapeChannel* lBlendShapeChannel =
FbxBlendShapeChannel::Create(pScene, "");

                FbxAnimCurve* myCompCurve =
lBlendShapeChannel->DeformPercent.GetCurve(myAnimBaseLayer, true);
                myCompCurve->KeyModifyBegin();
                FbxTime myTime;
                myTime.SetSecondDouble(0.0);
                int myKeyIndex = -1;

                myTime.SetSecondDouble((myKeyIndex + 1)*1.0);
                myKeyIndex = myCompCurve->KeyAdd(myTime);
                myCompCurve->KeySetValue(myKeyIndex, 0);
                myCompCurve->KeySetInterpolation(myKeyIndex,
FbxAnimCurveDef::eInterpolationLinear);

                std::string blendName("shape");
                blendName += std::to_string(i);

                FbxShape* shape = FbxShape::Create(pScene, blendName.c_str());
                shape->InitControlPoints(jointPos.cols());
                FbxVector4* pts = shape->GetControlPoints();
                DMatrix3X& inpts = global->componentRegressors[i];
                pts[0] = FbxVector4(0, 0, 0);
                for (int j = 0; j < global->joint_count; j++) {
                        pts[1 + j] = FbxVector4(
                                SCALAR * global->jointRegressor(3 * j + 0, 1 +
i).getValue(), SCALAR * global->jointRegressor(3 * j + 1, 1 + i).getValue(),
                                SCALAR * global->jointRegressor(3 * j + 2, 1 +
i).getValue()
                        );
                }

                for (int s = 0; s < global->persons.size(); s++) {
                        double coeff = global->persons[s]->blendCoefficient(0,
i).getValue(); for (int t = 0; t < global->persons[s]->scans.size(); t++) {

                                myTime.SetSecondDouble((myKeyIndex + 1)*1.0);
                                myKeyIndex = myCompCurve->KeyAdd(myTime);

                                myCompCurve->KeySetValue(myKeyIndex, 100 * coeff
/ SCALAR);

                                myCompCurve->KeySetInterpolation(myKeyIndex,
FbxAnimCurveDef::eInterpolationLinear);
                        }
                }

                lMesh->AddDeformer(lBlendShape);
                lBlendShape->AddBlendShapeChannel(lBlendShapeChannel);
                lBlendShapeChannel->AddTargetShape(shape);
                lBlendShapeChannel->DeformPercent.Set(0);
        }



        lNode->SetNodeAttribute(lMesh);
        return lNode;
}
*/
static FbxNode *createBaseMesh(FbxScene *pScene, const char *pName,
                               Model *model, int personId) {

  FbxMesh *lMesh = FbxMesh::Create(pScene, pName);

  FbxNode *lNode = FbxNode::Create(pScene, pName);

  // std::cout << "max: " << poseRegressor.maxCoeff() << " min: " <<
  // poseRegressor.minCoeff() << std::endl;

  Matrix3X a;
  a = model->mean;
  Person &person = model->persons[personId];
  // for (int j = 0; j < model->pcs.size(); j++) {
  //	a += person.coeff[j] * model->pcs[j];
  //}

  Eigen::Matrix<int, 5, Eigen::Dynamic> f = model->top.quads;
  lMesh->InitControlPoints(a.cols());
  FbxVector4 *lControlPoints = lMesh->GetControlPoints();

  for (int i = 0; i < model->mean.cols(); i++) {
    lControlPoints[i] = FbxVector4(a(0, i), a(1, i), a(2, i));
  }

  for (int i = 0; i < f.cols(); i++) {
    lMesh->BeginPolygon(-1, -1, -1, false);

    for (int j = 0; j < 4; j++) {
      lMesh->AddPolygon(f(j, i));
    }

    lMesh->EndPolygon();
  }
  /*lMesh->SetMeshSmoothness(FbxMesh::ESmoothness::eFine);
  lMesh->SetMeshPreviewDivisionLevels(3);
  lMesh->SetMeshRenderDivisionLevels(3);
  lMesh->SetDisplaySubdivisions(true);*/

  for (int i = 0; i < model->pcs.size(); i++) {
    for (int dir = -1; dir <= 1; dir += 2) {
      std::stringstream buffer;
      buffer << "shape" << std::setfill('0') << std::setw(5) << i;
      std::string blendName = buffer.str();
      FbxShape *shape = FbxShape::Create(pScene, blendName.c_str());
      shape->InitControlPoints(a.cols());
      FbxVector4 *pts = shape->GetControlPoints();
      int sign = (person.coeff[i] >= 0 ? 1 : -1);

      for (int j = 0; j < a.cols(); j++) {
        Vector3 pt = 10 * model->pcs[i].col(j);
        pts[j] = FbxVector4(dir * pt[0] + a(0, j), dir * pt[1] + a(1, j),
                            dir * pt[2] + a(2, j));
      }

      FbxBlendShape *lBlendShape = FbxBlendShape::Create(pScene, "");
      FbxBlendShapeChannel *lBlendShapeChannel =
          FbxBlendShapeChannel::Create(pScene, "");
      lMesh->AddDeformer(lBlendShape);
      lBlendShape->AddBlendShapeChannel(lBlendShapeChannel);
      lBlendShapeChannel->AddTargetShape(shape);
      // lBlendShapeChannel->DeformPercent.Set(10*fabs(person.coeff[i]));
      float mbeta = sign == dir ? 10 * fabs(person.coeff[i]) : 0;
      lBlendShapeChannel->DeformPercent.Set(mbeta);
      // lBlendShapeChannel->DeformPercent.Set(0);
    }
  }
  for (int i = 9; i < 9 * (model->joints.cols()) && model->deform.size() != 0;
       i++) {
    int jointIdx = i / 9;
    std::stringstream buffer;
    buffer << "pose" << std::setfill('0') << std::setw(5) << i;
    std::string blendName = buffer.str();

    FbxShape *shape = FbxShape::Create(pScene, blendName.c_str());
    shape->InitControlPoints(a.cols());
    FbxVector4 *pts = shape->GetControlPoints();
    // std::cout << "Model deform size: " << model->deform.size() << std::endl;
    // Matrix3X& inpts = model->deform[i];
    for (int j = 0; j < a.cols(); j++) {
      SVector4 influence = model->closestDeformJoints.col(j).cast<int>();
      // std::cout <<  "before find" << std::endl;

      // int* regIdx_ptr = std::find(&influence[0], &influence[0]+4, jointIdx);
      int *regIdx_ptr = std::find(&influence[0], &influence[0] + 2, jointIdx);

      // std::cout <<  "after find" << std::endl;
      int offset = regIdx_ptr - &influence[0];
      // std::cout <<  "offset: " << offset << std::endl;

      // std::cout << "offset: " << offset << std::endl;

      if (regIdx_ptr == &influence[0] + 2) {
        pts[j] = FbxVector4(a(0, j), a(1, j), a(2, j));
        // pts[j] = FbxVector4(0, 0, 0);
      } else {
        Vector3 pt = model->deform[9 * offset + (i % 9)].col(j);
        pts[j] = FbxVector4(pt[0] + a(0, j), pt[1] + a(1, j), pt[2] + a(2, j));
        if (j == 0) {

          // if (i % 9 == 0) {
          //	std::cout << blendName << " : " <<
          //model->skeleton->getJoint(jointIdx).second << 		" -> " <<
          //model->skeleton->getJoint(jointIdx).first << " : \n";
          //}
          // std::cout << pts[j][0] << ", "<< pts[j][1] << ", " << pts[j][2] <<
          // std::endl;
        }
      }
    }

    // std::cout << i << " : "  << pts[0][0] << ", " << pts[0][1] << ", " <<
    // pts[0][2] << std::endl;
    FbxBlendShape *lBlendShape = FbxBlendShape::Create(pScene, "");
    FbxBlendShapeChannel *lBlendShapeChannel =
        FbxBlendShapeChannel::Create(pScene, "");
    lMesh->AddDeformer(lBlendShape);
    lBlendShape->AddBlendShapeChannel(lBlendShapeChannel);
    lBlendShapeChannel->AddTargetShape(shape);
    lBlendShapeChannel->DeformPercent.Set(0);
  }
  lNode->SetNodeAttribute(lMesh);
  return lNode;
}
/*
static void dumpSubiv(std::string name, SubdivEvaluator* eval, const Matrix3X&
vertices) { MeshTopology resultTopology; Matrix3X resultVertices;
        eval->generate_refined_mesh(vertices, 1, &resultTopology,
&resultVertices); auto path = name + std::string(".obj"); saveObj(path,
&resultTopology, &resultVertices);
}

static void dumpMesh(std::string name, MeshTopology& topology, const Matrix3X&
vertices) { auto path = name + std::string(".obj"); saveObj(path, &topology,
&vertices);
}

static void dumpCloud(std::string name, int iter, const Matrix3X& cloud, const
Matrix3X* normals = NULL) { std::ofstream file(name + std::string(".xyz")); if
(normals == NULL) { file << cloud.transpose().format(XYZFormat);
        }
        else {
                MatrixX M(6, cloud.cols());
                M.block(0, 0, 3, cloud.cols()) = cloud;
                M.block(3, 0, 3, cloud.cols()) = *normals;
                file << M.transpose().format(XYZFormat);
        }
}
*/

bool reparameterize(const Model &source, Model *target) {
  MatrixX shapes;
  shapes.resize(source.scans.size(),
                3 * source.mean.cols() + 3 * source.joints.cols());

  for (int i = 0; i < source.mean.cols(); i++) {
    /*
    for(int j = 0; j < source.scans.size(); j++) {
        shapes(j,3*i+0) =  source.mean(0,i);
        shapes(j,3*i+1) =  source.mean(1,i);
        shapes(j,3*i+2) =  source.mean(2,i);
    }
    */

    for (int j = 0; j < source.scans.size(); j++) {
      for (int k = 0; k < source.pcs.size(); k++) {
        shapes(j, 3 * i + 0) += source.pcs[k](0, i) * 1;
        shapes(j, 3 * i + 1) += source.pcs[k](1, i) * 1;
        shapes(j, 3 * i + 2) += source.pcs[k](2, i) * 1;
      }
    }
  }
}
