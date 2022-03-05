#pragma once
#include <iostream>
#include "../CurveR2.hpp"
#include <queue>
#include "../DiscreteFrechet.hpp"


// Tree Node Structure.
class TreeNode{
    public:
        //Store Cuvre
        R2_Curve* curve;
        //Store Children
        TreeNode* left;
        TreeNode* right;
};

class Tree{
    public:
        TreeNode* root;


        //Generate Tree.
        Tree( vector<R2_Curve*> Curves){
            int number_of_curves = Curves.size();
            cout<<Curves.at(0)->get_name()<<endl;
            cout<<"Num_of_curves On tree"<<number_of_curves;
            queue<pair<int, TreeNode *>> Container; 
            
            while( number_of_curves > 0){

                // Generate Curve - Leaf - Node .//
                TreeNode* Curve_Leaf = new TreeNode;
                TreeNode* Left = NULL;
                TreeNode* Right = NULL;
                Curve_Leaf ->curve = Curves[ number_of_curves -1];
                // Keep The Nodes in a Container for Later - Use. //
                Container.emplace(1, Curve_Leaf);
                number_of_curves--;

            }

             cout<<"Container Size "<<Container.size()<<endl;

             while(1){
            
                //Get Nodes One by One And Generate Inner Nodes. //
                pair<unsigned int, TreeNode *>  CurveNode_1 = Container.front();
                Container.pop();
                cout<<"size"<<Container.size()<<endl;
                if (Container.size() == 0) 
                {       cout<<"edw"<<endl;
                        cout<<"Root "<<CurveNode_1.second->curve->get_name()<<endl;
                        root = CurveNode_1.second; 
                        break;
                }
                if (CurveNode_1.first < Container.front().first) 
                {
                     Container.push(CurveNode_1);
                }           
                else
                {
                    pair<int, TreeNode*> CurveNode_2 = Container.front();
                    Container.pop();

                    // Generate Inner Node Of 1-2 .//
                    TreeNode * InnerNode = new TreeNode;
                    InnerNode->left = CurveNode_1.second ;
                    InnerNode->right = CurveNode_2.second;
                    InnerNode->curve = NULL;
                    Container.emplace(CurveNode_1.first + CurveNode_2.first, InnerNode);
                }
            }


         }

    



        // Function to Print Tree .
        void PrintTree(TreeNode*  leaf)
        {
        if (leaf != NULL)
        {
            PrintTree(leaf->left);
            PrintTree(leaf->right);
            if ( leaf->curve != NULL)
            {
                cout<<leaf->curve->get_name()<<endl;
            }
        }

        }
















   /*

    // Post Order Traversal Implementation . 
    R2_Curve* PostOrderTraversal(TreeNode* node, vector<R2_Curve*> &temp_centroids){ 
        //If node is a Leaf . // 
        if (node->curve != NULL){
            return node->curve;
        }else{
            R2_Curve* LeftCurve = PostOrderTraversal(node->left, temp_centroids);
            R2_Curve* RightCurve;
            if (node->right != NULL){
                RightCurve = PostOrderTraversal(node->right, temp_centroids);
            }
            else{
                RightCurve = NULL;
            }

            R2_Curve* temp_centroid = MeanCurve(LeftCurve, RightCurve);
            node->curve = temp_centroid;
            if (RightCurve != NULL)
            {
                temp_centroids.push_back(temp_centroid); //keep temporally created curves into a vector so that they can been deleted later avoiding memory leaks
            }
            return temp_centroid;
        }
    }














    // Calculate Mean between 2 points .//
    Point CalculateMeanPoint(Point &point1, Point &point2){
        Point MeanPoint((point1.x_i + point2.x_i / 2.0) , (point2.y_i + point2.y_i / 2.0) );
        return MeanPoint;
    }


    void Get_Distances(const std::vector<Point> & v1,const std::vector<Point> &v2 , double **  &Distances){

        int n = v1.size();
        int m = v2.size();
        cout<<" curve_size ["<<n<<","<<"]"<<endl;
        double ** Calculation;
        Calculation = (double**) malloc(n * sizeof(double*));
        if (Calculation == NULL){
            printf("Malloc: memory allocation error!\n");
            exit(3);
        }
        for (int i = 0 ; i < n ; i++){
            Calculation[i] = (double*) malloc(m * sizeof(double));	
            if (Calculation[i] == NULL){
                printf("Malloc: memory allocation error!\n");
                exit(3);
            }
        }
        Calculation[0][0] = Euclidean_Distance(v1.at(0), v2.at(0));
        cout<<"Calculation[0][0] : "<<Calculation[0][0] <<endl;

        for (int i = 1 ; i < m ; i++)
            Calculation[i][0] = maximum(Calculation[i-1][0],  Euclidean_Distance(v1.at(i), v2.at(0)));
        
        for (int j = 1 ; j < n ; j++)
            Calculation[0][j] = maximum(Calculation[0][j-1],  Euclidean_Distance(v1.at(0), v2.at(j)));

        for (int i = 1 ; i < m ; i++)
            for (int j = 1 ; j < n ; j++)
                Calculation[i][j] = maximum(minimum_of_3(Calculation[i-1][j], Calculation[i][j-1], Calculation[i-1][j-1]), Euclidean_Distance(v1.at(i), v2.at(j)));


        Distances = Calculation;


        return;
    }











        
    R2_Curve * MeanCurve(R2_Curve * Curve_1, R2_Curve*  Curve_2){ 
        if ( Curve_2 == NULL)
            return Curve_1;
        if ( Curve_1 == NULL)
            return Curve_2;
    
        // Get Number of Points Per Curve .//
        unsigned int n = Curve_1->get_number_of_points();
        unsigned int m = Curve_2->get_number_of_points();

        //Get Their Points . //
        vector<Point> points_1 = Curve_1 ->getVector();
        vector<Point> points_2 = Curve_2 ->getVector();

        // Get Distances .//
        double ** Distances;
        Get_Distances(points_1,points_2, Distances);
        // Container to keep The mean points of the 2 Curves. //
        vector<Point> MeanPoints;
        MeanPoints.push_back(CalculateMeanPoint(points_1[n-1],points_2[m-1]));

        unsigned int i = n - 1 ;
        unsigned int j = m - 1;
        unsigned int index ;
        unsigned int index2;
        double temp_dfd ;

        while (1){
            if  ((i == 0) &&( j == 0)){
                break;
            }   
            if ( i == 0){
                j--;
            }
            else if ( j == 0){
                i--;
            }
            else{
                if( Distances[i-1][j] < Distances[i][j-1] )
                {
                    index2 = 0;
                    temp_dfd = Distances[i-1][j] ;
                }
            else
            {
                index2 = 1;
                temp_dfd = Distances[i][j-1] ;
            }

            index = temp_dfd < Distances[i-1][j-1] ? index2 : 2 ;

            if (index == 0)
            {
                i--;
            }
            else if (index == 1)
            {
                j--;
            }
            else if ( index == 2)
            {
                i--;
                j--;
            }
        }

        MeanPoints.push_back(CalculateMeanPoint(points_1[n-1],points_2[m-1]));
    }

    //vector<double> vector_double_null;
    R2_Curve* Mean_Curve =  new R2_Curve();
    for(int i = 0  ; i< MeanPoints.size() ; i++){
        Mean_Curve->setPoint(MeanPoints.at(i));
    }
    return Mean_Curve; 
}




    
    R2_Curve * GetMeanCurve(){
        vector<R2_Curve *> temp_centroids;
        R2_Curve*  NewCentroid =  PostOrderTraversal(root, temp_centroids);
        for(unsigned int i = 0; i < temp_centroids.size() - 1; i++){ 
            delete temp_centroids[i];
            temp_centroids[i] = NULL;
        }
        temp_centroids.clear();

        return NewCentroid;
}   





 





*/





};