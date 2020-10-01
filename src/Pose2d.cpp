#include "Pose2d.h"

#include <fstream>
#include <nlohmann/json.hpp>
#include <map>
#include <iostream>

using json = nlohmann::json;
static const int H135=26;
static const int R135=27;
static const int F135 = R135+40;
static std::map<unsigned int, std::string> POSE_BODY_25_BODY_PARTS {
     {0,  "Nose"},
     {1,  "Neck"},
     {2,  "RShoulder"},
     {3,  "RElbow"},
     {4,  "RWrist"},
     {5,  "LShoulder"},
     {6,  "LElbow"},
     {7,  "LWrist"},
     {8,  "MidHip"},
     {9,  "RHip"},
     {10, "RKnee"},
     {11, "RAnkle"},
     {12, "LHip"},
     {13, "LKnee"},
     {14, "LAnkle"},
     {15, "REye"},
     {16, "LEye"},
     {17, "REar"},
     {18, "LEar"},
     {19, "LBigToe"},
     {20, "LSmallToe"},
     {21, "LHeel"},
     {22, "RBigToe"},
     {23, "RSmallToe"},
     {24, "RHeel"},

	//Left Hand
	{25,"LWrist2"},
	{H135+0, "LThumb1CMC"},       {H135+1, "LThumb2Knuckles"}, {H135+2, "LThumb3IP"},   {H135+3, "LThumb4FingerTip"},
    {H135+4, "LIndex1Knuckles"},  {H135+5, "LIndex2PIP"},      {H135+6, "LIndex3DIP"},  {H135+7, "LIndex4FingerTip"},
    {H135+8, "LMiddle1Knuckles"}, {H135+9, "LMiddle2PIP"},     {H135+10, "LMiddle3DIP"},{H135+11, "LMiddle4FingerTip"},
    {H135+12, "LRing1Knuckles"},  {H135+13, "LRing2PIP"},      {H135+14, "LRing3DIP"},  {H135+15, "LRing4FingerTip"},
	{H135+16, "LPinky1Knuckles"}, {H135+17, "LPinky2PIP"}, {H135+18, "LPinky3DIP"}, {H135+19, "LPinky4FingerTip"},

	//Right Hand
	{46,"RWrist2"},
	{R135+20, "RThumb1CMC"},      {R135+21, "RThumb2Knuckles"},{R135+22, "RThumb3IP"},  {R135+23, "RThumb4FingerTip"},
    {R135+24, "RIndex1Knuckles"}, {R135+25, "RIndex2PIP"},     {R135+26, "RIndex3DIP"}, {R135+27, "RIndex4FingerTip"},
    {R135+28, "RMiddle1Knuckles"},{R135+29, "RMiddle2PIP"},    {R135+30, "RMiddle3DIP"},{R135+31, "RMiddle4FingerTip"},
    {R135+32, "RRing1Knuckles"},  {R135+33, "RRing2PIP"},      {R135+34, "RRing3DIP"},  {R135+35, "RRing4FingerTip"},
   	{R135+36, "RPinky1Knuckles"}, {R135+37, "RPinky2PIP"}, {R135+38, "RPinky3DIP"}, {R135+39, "RPinky4FingerTip"},

	{F135+0, "FaceContour0"},   {F135+1, "FaceContour1"},   {F135+2, "FaceContour2"},   {F135+3, "FaceContour3"},   {F135+4, "FaceContour4"},   {F135+5, "FaceContour5"},   // Contour 1/3
    {F135+6, "FaceContour6"},   {F135+7, "FaceContour7"},   {F135+8, "FaceContour8"},   {F135+9, "FaceContour9"},   {F135+10, "FaceContour10"}, {F135+11, "FaceContour11"}, // Contour 2/3
    {F135+12, "FaceContour12"}, {F135+13, "FaceContour13"}, {F135+14, "FaceContour14"}, {F135+15, "FaceContour15"}, {F135+16, "FaceContour16"},                             // Contour 3/3
    {F135+17, "REyeBrow0"},  {F135+18, "REyeBrow1"},  {F135+19, "REyeBrow2"},  {F135+20, "REyeBrow3"},  {F135+21, "REyeBrow4"}, // Right eyebrow
    {F135+22, "LEyeBrow4"},  {F135+23, "LEyeBrow3"},  {F135+24, "LEyeBrow2"},  {F135+25, "LEyeBrow1"},  {F135+26, "LEyeBrow0"}, // Left eyebrow
    {F135+27, "NoseUpper0"}, {F135+28, "NoseUpper1"}, {F135+29, "NoseUpper2"}, {F135+30, "NoseUpper3"}, // Upper nose
    {F135+31, "NoseLower0"}, {F135+32, "NoseLower1"}, {F135+33, "NoseLower2"}, {F135+34, "NoseLower3"}, {F135+35, "NoseLower4"}, // Lower nose
    {F135+36, "REye0"}, {F135+37, "REye1"}, {F135+38, "REye2"}, {F135+39, "REye3"}, {F135+40, "REye4"}, {F135+41, "REye5"}, // Right eye
    {F135+42, "LEye0"}, {F135+43, "LEye1"}, {F135+44, "LEye2"}, {F135+45, "LEye3"}, {F135+46, "LEye4"}, {F135+47, "LEye5"}, // Left eye
    {F135+48, "OMouth0"}, {F135+49, "OMouth1"}, {F135+50, "OMouth2"}, {F135+51, "OMouth3"}, {F135+52, "OMouth4"}, {F135+53, "OMouth5"}, // Outer mouth 1/2
    {F135+54, "OMouth6"}, {F135+55, "OMouth7"}, {F135+56, "OMouth8"}, {F135+57, "OMouth9"}, {F135+58, "OMouth10"}, {F135+59, "OMouth11"}, // Outer mouth 2/2
    {F135+60, "IMouth0"}, {F135+61, "IMouth1"}, {F135+62, "IMouth2"}, {F135+63, "IMouth3"}, {F135+64, "IMouth4"}, {F135+65, "IMouth5"}, {F135+66, "IMouth6"}, {F135+67, "IMouth7"}, // Inner mouth
	{F135+68, "RPupil"}, {F135+69, "LPupil"} // Pupils
};

static std::map<unsigned int, std::string> POSE_BODY_18_BODY_PARTS {
     {0,  "Nose"},
     {1,  "Neck"},
     {2,  "RShoulder"},
     {3,  "RElbow"},
     {4,  "RWrist"},
     {5,  "LShoulder"},
     {6,  "LElbow"},
     {7,  "LWrist"},
     {8,  "RHip"},
     {9,  "RKnee"},
     {10, "RAnkle"},
     {11, "LHip"},
     {12, "LKnee"},
     {13, "LAnkle"},
     {14, "REye"},
     {15, "LEye"},
     {16, "REar"},
     {17, "LEar"},
};


static const std::vector<std::string> part_names_ar = {
                "m_chest",
                "m_neck",
                "m_head",
                "m_rupperarm",
                "m_rlowerarm",
                "m_rhand",
                "m_rthumb",
                "m_rthumb-proximal",
                "m_rthumb-distal",
                "m_rfist",
                "m_rindex-middle",
                "m_rmiddle-middle",
                "m_rring-middle",
                "m_rpinky-middle",
                "m_lupperarm",
                "m_llowerarm",
                "m_lhand",
                "m_lthumb",
                "m_lthumb-proximal",
                "m_lthumb-distal",
                "m_lfist",
                "m_lindex-middle",
                "m_lmiddle-middle",
                "m_lring-middle",
                "m_lpinky-middle",
                "m_hip",
                "m_pelvis",
                "m_rupperleg",
                "m_rlowerleg",
                "m_rfoot",
                "m_lupperleg",
                "m_llowerleg",
                "m_lfoot"
            };




static const std::map<std::string, std::string> openpose_to_skel {
/*    {"m_neck","Neck"},
    {"m_rupperarm","LShoulder"},
    {"m_rlowerarm","LElbow"},
    {"m_rhand","LWrist"},
    {"m_lupperarm","RShoulder"},
    {"m_llowerarm","RElbow"},
    {"m_lhand","RWrist"},
    {"m_rupperleg","LHip"},
    {"m_rlowerleg","LKnee"},
    {"m_rfoot","LAnkle"},
    {"m_lupperleg","RHip"},
    {"m_llowerleg","RKnee"},
    {"m_lfoot","RAnkle"},

    {"m_rthumb-proximal","LThumb2Knuckles"},
    {"m_rthumb-distal","LThumb3IP"},
    {"m_rindex-middle","LIndex2PIP"},
    {"m_rmiddle-middle","LMiddle2PIP"},
    {"m_rring-middle","LRing2PIP"},
    {"m_rpinky-middle","LPinky2PIP"},

	{"m_lthumb-proximal","RThumb2Knuckles"},
    {"m_lthumb-distal","RThumb3IP"},
    {"m_lindex-middle","RIndex2PIP"},
    {"m_lmiddle-middle","RMiddle2PIP"},
    {"m_lring-middle","RRing2PIP"},
    {"m_lpinky-middle","RPinky2PIP"}
*/
};

static std::map<size_t, std::string> surface_to_openpose {
/*
    {1, "Nose"},
	{658, "RBigToe"},
	{637, "RSmallToe"},
	{222, "LBigToe"},
	{198, "LSmallToe"},
	{1141,"RThumb4FingerTip"},
	{876,"RIndex4FingerTip"},
	{904,"RMiddle4FingerTip"},
	{874,"RRing4FingerTip"},
	{907,"RPinky4FingerTip"},

    {1165,"LThumb4FingerTip"},
    {1011,"LIndex4FingerTip"},
    {1038,"LMiddle4FingerTip"},
    {1320,"LRing4FingerTip"},
    {1041,"LPinky4FingerTip"},

	{522, "REyeBrow0"},  {521, "REyeBrow1"},  {506, "REyeBrow2"},  {505, "REyeBrow3"},
	{66, "LEyeBrow3"},  {67, "LEyeBrow2"},  {82, "LEyeBrow1"},  {83, "LEyeBrow0"},


	{7, "NoseUpper0"}, {4, "NoseUpper1"}, {3, "NoseUpper2"}, 
	{29, "NoseLower4"}, {30, "NoseLower2"}, {2, "NoseLower2"}, {469, "NoseLower1"}, {468, "NoseLower0"},


	{568, "OMouth0"}, {566, "OMouth2"}, {134, "OMouth3"}, {127, "OMouth4"}, 
    {129, "OMouth6"}, {132, "OMouth8"}, {133, "OMouth9"},{571, "OMouth10"}

*/
};


std::string Pose2d::convertJointName(std::string joint) {
	auto it = openpose_to_skel.find(joint);
	if(it != openpose_to_skel.end())
	{
   		return it->second;
	}
	return "";
}

Pose2d::Pose2d(std::string surfaceMapPath) {
	if(surface_to_openpose.size() > 0)
		return;

	std::ifstream file(surfaceMapPath);
    json j;
    file >> j;
	for (json::iterator it = j.begin(); it != j.end(); ++it) {
    	std::cout << it.key() << " has value " << it.value() << std::endl;
		surface_to_openpose[std::stoi(it.key())] = it.value();
	}
}

void Pose2d::init(bool base25) {
	auto& POSE = base25 ? POSE_BODY_25_BODY_PARTS : POSE_BODY_18_BODY_PARTS; 
	joints.resize(3,POSE.size());
	joints.setConstant(0);

	{
    	for(auto forwardIT = openpose_to_skel.begin();
            forwardIT != openpose_to_skel.end();
            forwardIT++)
    	{
    	    skel_to_openpose[forwardIT->second] = forwardIT->first; 
    	}
	}

	{
    	for(auto forwardIT = POSE.begin(); 
    	        forwardIT != POSE.end(); 
    	        forwardIT++)
   		{
	       	openpose_to_index[forwardIT->second] = forwardIT->first;
    	}
	}

}

Vector3 Pose2d::getJoint(std::string joint) {
	if(joints.cols() == 0) return Vector3(0,0,0);
	std::string openposeJoint = convertJointName(joint);
	if(openposeJoint == "") {
		Vector3 j; j.setConstant(0);
		return j;
	}
	unsigned int index = openpose_to_index[openposeJoint];
	std::cout << "Joint2d Index: " << index << std::endl;
	Vector3 result =  joints.col(index);
	//if(result(2) < .0)
    //    result(2) = 0;
	return result;
}

std::vector<size_t> Pose2d::getSurfaceIndices() {
	std::vector<size_t> keys;
	std::transform(
    surface_to_openpose.begin(),
    surface_to_openpose.end(),
    std::back_inserter(keys),
    [](const std::map<size_t,std::string>::value_type &pair){return pair.first;});
	return keys;
}

Vector3 Pose2d::getSurface(size_t vertexIdx) {
	if(joints.cols() == 0) return Vector3(0,0,0);
	if ( surface_to_openpose.find(vertexIdx) == surface_to_openpose.end() )
		return Vector3(0,0,0);
	std::string name = surface_to_openpose[vertexIdx];
	size_t openpose_idx = openpose_to_index[name];
	Vector3 result = joints.col(openpose_idx);
	//if(result(2) < .0)
    //    result(2) = 0;

	if(openpose_idx >= F135) {
		result(2) *= 5.0;
	}
	return result;
}

void Pose2d::read(std::string path) {
	std::cout << "Openpose Version: " << path << std::endl;
	std::ifstream file(path);
	json j;
	file >> j;
	std::cout << j["version"].get<float>()	 << std::endl;
	if(j["version"].get<float>() >= 1.2f) {
		init(true);
		if(j["people"].size() == 0) return;
		auto raw_joints =  j["people"][0]["pose_keypoints_2d"];
		auto left_hand = j["people"][0]["hand_left_keypoints_2d"];
        auto right_hand = j["people"][0]["hand_right_keypoints_2d"];
		auto head = j["people"][0]["face_keypoints_2d"];
        raw_joints.insert(raw_joints.end(), left_hand.begin(), left_hand.end());
        raw_joints.insert(raw_joints.end(), right_hand.begin(), right_hand.end());
		raw_joints.insert(raw_joints.end(), head.begin(), head.end());
		std::vector<Scalar> vec_joints;
		raw_joints.get_to(vec_joints);
		assert((POSE_BODY_25_BODY_PARTS.size()) * 3 == vec_joints.size());
		for(int i = 0; i < vec_joints.size(); i++) {
			joints(i%3,i/3) = vec_joints[i];
		}
	} else if(j["version"].get<float>() == 1.0f) {
		assert(false);
		init(false);
		if(j["people"].size() == 0) return;
		auto raw_joints =  j["people"][0]["pose_keypoints_2d"];
	
        std::vector<Scalar> vec_joints;
        
		raw_joints.get_to(vec_joints);
        assert((POSE_BODY_18_BODY_PARTS.size()) * 3 == vec_joints.size());
		
		
        for(int i = 0; i < vec_joints.size(); i++) {
            joints(i%3,i/3) = vec_joints[i];
        }


	} else {
		assert(false);
	}
}
