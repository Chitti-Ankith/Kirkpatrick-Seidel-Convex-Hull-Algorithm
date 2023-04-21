#include<bits/stdc++.h>
#include "Point.h"

using namespace std;

void swapd(double* a,double* b){
	double t;
	t = *a;
	*a=*b;
	*b=t;
}


bool leftu(double p1, double p2)
{
    return (p1 < p2);
}

bool leftpu(Point p1, Point p2)
{
    return (p1.getX() < p2.getX() || (p1.getX() == p2.getX() && p1.getY() > p2.getY()));
}

bool rightu(Point p1, Point p2)
{
    return (p1.getX() > p2.getX() || (p1.getX() == p2.getX() && p1.getY() > p2.getY()));
}

Point leftmostPointu(vector<Point> points)
{
    sort(points.begin(), points.end(), leftpu);
    return points[0];
}

Point rightmostPointu(vector<Point> points)
{
    sort(points.begin(), points.end(), rightu);
    return points[0];
}

int partitiond(vector<double> &points,int l,int r,double x)
{
    int i;
    for (i=l; i<r; i++)
        if (points[i] == x)
           break;
    swapd(&points[i],&points[r]);

    i = l;
    for (int j = l; j <= r - 1; j++)
    {
        if (points[j] <= x)
        {
            swapd(&points[j], &points[i]);
            i++;
        }
    }
    swapd(&points[i], &points[r]);
    return i;
}

double momd(vector<double> &points,int l,int r,int k)
{

        int num_points = r-l+1; // Number of elements in arr[l..r]

        // Divide arr[] in groups of size 5, calculate median
        // of every group and store it in median[] array.
        //int i;
        vector<double> median; // There will be floor((n+4)/5) groups;
        /*for (i=l; i<=r; i += 5)
            median[i] = findMed(points[l+i*5], 5);
        if (i*5 < n) //For last group with less than 5 elements
        {
            median[i] = findMed(points[l+i*5], n%5);
            i++;
        }*/

        for(int i = l;i<=r;i+=5)
        {
          int start=i;
          int end;
          if(i+5 <=r)
            end = i+5;
          else
            end = r+1;

      		sort(points.begin()+start,points.begin()+end,leftu);
          median.push_back(points[(start+end)/2]);
        }
        // Find median of all medians using recursive call.
        // If median[] has only one element, then no need
        // of recursive call
        double medianofmedian = (median.size()==1) ? median[0]:
                                 momd(median, 0, median.size()-1, median.size()/2);

        // Partition the array around a random element and
        // get position of pivot element in sorted array
        int pos = partitiond(points, l, r, medianofmedian);

        // If position is same as k
        if (pos-l == k)
            return points[pos];
        if (pos-l > k)  // If position is more, recur for left
            return momd(points, l, pos-1, k);

        // Else recur for right subarray
        return momd(points, pos+1, r, k-pos+l-1);

}

pair<Point,Point> upperBridge(vector<Point> S,int num_points,Point L,int depth)
{
  cout << "Upper Bridge Depth" << depth << endl;
  // cout << S.size() << endl;
  // cout << num_points << endl;
	if(S.size()==2)
  {
    if(S[0].getX() <= S[1].getX())
		return make_pair(S[0],S[1]);
    else
    return make_pair(S[1],S[0]);
  }
  vector<Point> candidates;
  vector<pair<Point,Point>> pairs;

  random_shuffle(S.begin(), S.end());
	
	for(int i =0;i< S.size();i++)
	S[i].printPoint();

  if(num_points%2==1)
  {
	{
		candidates.push_back(S[num_points-1]);
	}
  
  for(int i=0;i<num_points-1;i=i+2)
	{
    if(S[i].getX() <= S[i+1].getX())
		pairs.push_back(make_pair(S[i],S[i+1]));
    else
    pairs.push_back(make_pair(S[i+1],S[i]));
	}
  }
  else
  {
    
	  for(int i=0;i<num_points;i=i+2)
	{
		if(S[i].getX() <= S[i+1].getX())
		pairs.push_back(make_pair(S[i],S[i+1]));
    else
    pairs.push_back(make_pair(S[i+1],S[i]));
	}
  }

	//map<pair<Point, Point>, double> K;
  vector<double> kmom; //Simplifies implementation of medianofmedians
  vector<double> km;
	for(int i=0;i<pairs.size();i++)
	{	
		
		if(pairs[i].first.getX()==pairs[i].second.getX())  //Slope is 0
		{
			if(pairs[i].first.getY()>pairs[i].second.getY())
				candidates.push_back(pairs[i].first);
			else
				candidates.push_back(pairs[i].second);
			//K.push_back(make_pair(0, make_pair(pairs[i].first, pairs[i].second)));
      //K.insert(make_pair(make_pair(pairs[i].first,pairs[i].second), 0));
      //kmom.push_back(0);
      km.push_back(10000);
		}
		else
		{
      double slope = (pairs[i].first.getY()-pairs[i].second.getY())/(pairs[i].first.getX()-pairs[i].second.getX());
			//K.insert(make_pair(make_pair(pairs[i].first,pairs[i].second), slope));
      kmom.push_back(slope);
      km.push_back(slope);
			//cout<<K[i]<<endl;
		}
	}

  vector<double> kmomc;
  kmomc = kmom;
	//  for(int i =0;i< kmom.size();i++)
	//  {
	// 	 cout << kmom[i] << endl;
     
	// 	//  pairs[i].first.printPoint();
	// 	// // pairs[i].second.printPoint();
		
	//  }
  //  for(int i=0;i<km.size();i++)
  //  cout << "Km" << km[i] << endl;
 // need to define another medianofmedians for the double values called momd(below). Should be a copy/paste
  double k = momd(kmomc,0,kmomc.size()-1,kmomc.size()/2); 
  vector<pair<Point,Point>> SMALL;
  vector<pair<Point,Point>> EQUAL;
  vector<pair<Point,Point>> LARGE;
  // cout<<"mean slope ="<<k<<endl;
  for(int i=0;i<pairs.size();i++)
  { 
    // cout <<i<< "  "<<kmom[i]<<endl;
    if(km[i] != 10000)
    {
    if(k < km[i])
    {
      //  cout<<i<<" LSARGE "<<k<<" "<<kmom[i]<<endl;
    LARGE.push_back(make_pair(pairs[i].first,pairs[i].second));
    }
    else if(k == km[i])
    {
      // cout<<i<<" equal "<<k<<" "<<kmom[i]<<endl;
    EQUAL.push_back(make_pair(pairs[i].first,pairs[i].second));
    
    }
    else
    {
      // cout<<i<<" small "<<k<<" "<<kmom[i]<<endl;
    SMALL.push_back(make_pair(pairs[i].first,pairs[i].second));
    }
    }
  }

  int max = -10000;
  vector<Point> MAX;
  int y_intersection[num_points];
  for(int i=0;i<num_points;i++)
  {
    y_intersection[i] = S[i].getY() - k*S[i].getX();
    if(y_intersection[i] > max)
    {
      max = y_intersection[i];

    }
  }
  // cout<<"max::"<<max<<endl;

  for(int i=0;i<num_points;i++)
  {
    if(y_intersection[i] == max)
    MAX.push_back(S[i]);
  }
  // cout<<"max size"<<MAX.size()<<endl;
  // for(int i=0;i<MAX.size();i++)
  // {cout<<i<<"max";
  //   MAX[i].printPoint();
  // }
  Point pk = leftmostPointu(MAX);
  Point pm = rightmostPointu(MAX);
 
  if(pk.getX() <= L.getX() && pm.getX() > L.getX())
  {
  cout<<"pk is :";
  pk.printPoint();
  cout<<"pm is :";
  pm.printPoint();
  cout<<endl;
  return make_pair(pk,pm);

  }
  if(pm.getX() <= L.getX())
  {
    for(int i=0;i<LARGE.size();i++)
    candidates.push_back(LARGE[i].second);

    for(int i=0;i<EQUAL.size();i++)
    candidates.push_back(EQUAL[i].second);

    for(int i=0;i<SMALL.size();i++)
    {
      candidates.push_back(SMALL[i].first);
      candidates.push_back(SMALL[i].second);
    }

  }

  if(pk.getX() > L.getX())
  {
    for(int i=0;i<LARGE.size();i++)
    {
      candidates.push_back(LARGE[i].first);
      candidates.push_back(LARGE[i].second);
    }
    for(int i=0;i<EQUAL.size();i++)
    candidates.push_back(EQUAL[i].first);

    for(int i=0;i<SMALL.size();i++)
    {
      candidates.push_back(SMALL[i].first);
    }

  }

  // cout << "Hello\n"<<candidates.size();
  // for(int i =0;i< candidates.size();i++)
	// candidates[i].printPoint();
  //sort(candidates.begin(),candidates.end(),leftpu);
  return upperBridge(candidates,candidates.size(),L,depth+1);

}

// ***************************************Upper HUll Part*******************************

void swap(Point* a,Point* b){
	Point t;
	t = *a;
	*a=*b;
	*b=t;
}


int partition(vector<Point> &points,int l,int r,Point x)
{
    int i;
    for (i=l; i<r; i++)
        if (points[i].getX() == x.getX() && points[i].getY() == x.getY())
           break;
    swap(&points[i],&points[r]);

    i = l;
    for (int j = l; j <= r - 1; j++)
    {
        if (points[j].getX() <= x.getX())
        {
            swap(&points[j], &points[i]);
            i++;
        }
    }
    swap(&points[i], &points[r]);
    return i;
}



Point mom(vector<Point> &points,int l,int r,int k)
{

        int num_points = r-l+1; // Number of elements in arr[l..r]

        vector<Point> median; 
        /*for (i=l; i<=r; i += 5)
            median[i] = findMed(points[l+i*5], 5);
        if (i*5 < n) //For last group with less than 5 elements
        {
            median[i] = findMed(points[l+i*5], n%5);
            i++;
        }*/

        for(int i = l;i<=r;i+=5)
        {
          int start=i;
          int end;
          if(i+5 <=r)
            end = i+5;
          else
            end = r+1;

      		sort(points.begin()+start,points.begin()+end,leftpu);
          median.push_back(points[(start+end)/2]);
        }
        // Find median of all medians using recursive call.
        // If median[] has only one element, then no need
        // of recursive call
        Point medianofmedian = (median.size()==1) ? median[0]:
                                 mom(median, 0, median.size()-1, median.size()/2);

        // Partition the array around a random element and
        // get position of pivot element in sorted array
        int pos = partition(points, l, r, medianofmedian);

        // If position is same as k
        if (pos-l == k)
            return points[pos];
        if (pos-l > k)  // If position is more, recur for left
            return mom(points, l, pos-1, k);

        // Else recur for right subarray
        return mom(points, pos+1, r, k-pos+l-1);

}

int directionOfPoint(Point A, Point B, Point C) 
{ 
  
    // // Determining cross Product 
    // int cross_product = (B.getX() - A.getX()) * (C.getY() - A.getY()) - (B.getY() - A.getY()) * (C.getX() - A.getX()); 
  
    // // return RIGHT if cross product is positive 
    // if (cross_product > 0) 
    //     return 1; 
  
    // // return LEFT if cross product is negative 
    // if (cross_product < 0) 
    //     return -1; 
  
    // // return ZERO if cross product is zero.  
    // return 0; 

    double slope = (A.getY()-B.getY())/(A.getX()-B.getX());
    double value = slope*C.getX();
    double y0 = value + A.getY()- value*A.getX();

  if(slope>0)
  {
    if(C.getY()-y0 > 0)
    {
      return -1;
    }
    else if(C.getY()-y0 < 0)
    {
      return 1;
    }
    else
    {
      return 0;
    }
  }
  else
  {
    if(C.getY()-y0 > 0)
    {
      cout<<"point to the right"<<C.getX()<<C.getY()<<endl;
      return 1;
    }
    else if(C.getY()-y0 < 0)
    {
      return -1;
    }
    else
    {
      return 0;
    }
  }
  

} 


vector<Point> Upperhull(Point pmin,Point pmax,vector<Point> points,int depth)
{

    cout << "Depth =" << depth;
    // pmin.printPoint();
    // pmax.printPoint();
    if(pmin.getX() == pmax.getX() && pmax.getX() == pmax.getY())
    {
      vector<Point> res;
      res.push_back(pmin);
      return res;
    }
    int num_points = points.size();
    Point a = mom(points,0,num_points-1,num_points/2);
    // cout << "Mom" << endl;
    // a.printPoint();

    vector<Point> TL,TR;
    for(int i=0;i<num_points;i++)
  	{
  		if(points[i].getX()<a.getX())
  		{
  			TL.push_back(points[i]);
  		}
  		else if(points[i].getX()>=a.getX())
  		{
  			TR.push_back(points[i]);
  		}
  	}

    // cout << "TL" << endl;

    // for(int i =0;i<TL.size();i++)
    // TL[i].printPoint();

    // cout << "TR" << endl;

    // for(int i =0;i<TR.size();i++)
    // TR[i].printPoint();

    pair<Point,Point> ub = upperBridge(points,num_points,a,1);
    Point pl,pr;
    pl = ub.first;
    pr = ub.second;

    cout << "Upper Bridge" << endl;
    pl.printPoint();
    pr.printPoint();

    vector<Point> Tl,Tr;
    Tl.push_back(pl);
    Tr.push_back(pr);
	
	//Need to add points to Tl and Tr as given in slides. Please Check This. Do we need to add points on that line as well?
	if(pl.getX() != pmin.getX() && pl.getY() != pmin.getY())
  {
    for(int i=0;i<TL.size();i++)
	{
		if(TL[i].getX() != pl.getX() || TL[i].getY() != pl.getY()) //since pl is already added
	    {
			if(directionOfPoint(pl,pmin,TL[i]) == -1)
			Tl.push_back(TL[i]);
		}
	}
  }
	for(int i=0;i<TR.size();i++)
	{
		if(TR[i].getX() != pr.getX() || TR[i].getY() != pr.getY())    
		{
			if(directionOfPoint(pr,pmax,TR[i]) == 1) //since pr is already added
			Tr.push_back(TR[i]);
		}
	}

  cout << "Tl" << endl;
  for(int i=0;i<Tl.size();i++)
    Tl[i].printPoint();

  cout << "Tr" << endl;
  for(int i=0;i<Tr.size();i++)
    Tr[i].printPoint();

  cout << "End of Tr" << endl;

    vector<Point> Uh1,Uh2;
    vector<Point> ansHull;
    cout << "Tl" << Tl.size() << endl;
    cout << "Tr" << Tr.size() << endl;
    if(Tl.size() >= 2)
    Uh1 = Upperhull(pmin,pl,Tl,depth+1);

    if(Tr.size() >= 2)
    Uh2 = Upperhull(pr,pmax,Tr,depth+1);

    ansHull.push_back(pmin);
    for(int i=0;i<Uh1.size();i++)
    ansHull.push_back(Uh1[i]);

    ansHull.push_back(pl);
    ansHull.push_back(pr);

    for(int i=0;i<Uh2.size();i++)
    ansHull.push_back(Uh2[i]);

    ansHull.push_back(pmax);
     //Concatenating the two vectors

    return ansHull;
}



int main(void)
{
    ifstream input("points.txt");
    double x, y;
    char comma;
    vector<Point> points;
    while (input >> x >> comma >> y)
    {
        points.push_back(Point(x,y));
    }
	//   Point L(34,67);
	// pair<Point,Point> ub = upperBridge(points,points.size(),L);
	// cout<<"point"<<endl;
	// ub.first.printPoint();
	// // cout<<"\n";
	// ub.second.printPoint();

    Point pumin,pumax;

    pumin = leftmostPointu(points);
    pumax = rightmostPointu(points);


    vector<Point> Tu;
  	vector<Point> Tl;
  	Tu.push_back(pumin);
  	Tu.push_back(pumax);

    for(int i=0;i<points.size();i++)
  	{
  		if(points[i].getX()>pumin.getX() && points[i].getX()<pumax.getX())
  		{
  			Tu.push_back(points[i]);
  		}
  	}



    vector<Point> uh = Upperhull(pumin,pumax,Tu,0);
    // cout << uh.size();
    // for(int i =0;i<uh.size();i++)
    // {
    //   uh[i].printPoint();
    // }

  //  sort( uh.begin(), uh.end(),leftpu );
  //   uh.erase( unique( uh.begin(), uh.end() ), uh.end() );

    for(int i=0;i<uh.size();i++)
    uh[i].printPoint();

}