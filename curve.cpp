#include "curve.h"
#include "extra.h"
#ifdef WIN32
#include <windows.h>
#endif
#include <GL/gl.h>
#include <fstream>
using namespace std;

 // Bernstein basis matrix
 const   Matrix4f bernstein(
                            1.0f, -3.0f, 3.0f, -1.0f,
                            0.0f, 3.0f, -6.0f, 3.0f,
                            0.0f, 0.0f, 3.0f, -3.0f,
                            0.0f, 0.0f, 0.0f, 1.0f
                        );

 // first derivative of bernstein basis matrix
 const   Matrix4f d_bernstein(
                            0.0f, -3.0f, 6.0f, -3.0f,
                            0.0f, 3.0f, -12.0f, 9.0f,
                            0.0f, 0.0f, 6.0f, -9.0f,
                            0.0f, 0.0f, 0.0f, 3.0f
                            );

 // b-spline basis matrix, all entries are divided by 6
 const Matrix4f b_spline(
    1.0/6, -3.0/6, 3.0/6, -1.0/6,
    4.0/6, 0.0/6, -6.0/6, 3.0/6,
    1.0/6, 3.0/6, 3.0/6, -3.0/6,
    0.0/6, 0.0/6, 0.0/6, 1.0/6
);

namespace
{
    // Approximately equal to.  We don't want to use == because of
    // precision issues with floating point.
    inline bool approx( const Vector3f& lhs, const Vector3f& rhs )
    {
        const float eps = 1e-8f;
        return ( lhs - rhs ).absSquared() < eps;
    }

    
}
    

Curve evalBezier( const vector< Vector3f >& P, unsigned steps )
{
    // Check
    if( P.size() < 4 || P.size() % 3 != 1 )
    {
        cerr << "evalBezier must be called with 3n+1 control points." << endl;
        exit( 0 );
    }

    Curve curve;

    float t;

    // initial seed for the binormals
    Vector3f B_prev(0.0f, 0.0f, 1.0f);

    // geometry matrix consisting of the coordinates of the 4 control points
    Matrix4f geo_matrix(
                            P[0].x(), P[1].x(), P[2].x(), P[3].x(),
                            P[0].y(), P[1].y(), P[2].y(), P[3].y(),
                            P[0].z(), P[1].z(), P[2].z(), P[3].z(),
                            0.0f, 0.0f, 0.0f, 1.0f
                            );

    

    cout << "\t>>> evalBezier has been called with the following input:" << endl;

    cout << "\t>>> Control points (type vector< Vector3f >): "<< endl;
    for( unsigned i = 0; i < P.size(); ++i )
    {
        cout << "\t>>> (" << P[i].x() << ", " << P[i].y() << ", " << P[i].z() << ")" << endl;
    }

    cout << "\t>>> Steps (type steps): " << steps << endl;

    // for each step in steps
    for(unsigned k = 1; k <= steps; k++) 
    {
        // the point we are plotting
        CurvePoint p;

        // derivative of the curve at t
        Vector4f d_q_t;

        // tangent vector
        Vector4f T;
        
        // t = i / steps
        t = static_cast<float>(k) / static_cast<float>(steps);

        // monomial basis: 1, t, t^2, t^3
        Vector4f monomials(1, t, pow(t, 2), pow(t, 3));
        
        //cout << "the value of t is: " << t << endl;

        //cout << "monomial basis: ";

        //monomials.print();

        // calculate length 4 vector q(t)
        Vector4f result = geo_matrix * bernstein * monomials;

        if(result.z() == -0)
            result.z() = 0.0;

        // store result coordinates in proper length 3 vector
        p.V = Vector3f(result.x(), result.y(), result.z());

        // calculate q'(t)
        d_q_t = geo_matrix * d_bernstein * monomials;

        //cout << "tangent before being normalized: ";

        //d_q_t.print();

        // Tangent = ||q'(t)||
        T = d_q_t.normalized();

        //cout << "tangent after being normalized: ";

        //T.print();

        // change any -0 z values to 0
        if(T.z() == -0)
            T.z() = 0.0;

        // store T
        p.T = Vector3f(T.x(), T.y(), T.z());

        // N_i = (B_i-1 * T_i).normalized
        p.N = (Vector3f::cross(B_prev, p.T)).normalized();

        // get rid of -0 values
        if(p.N.z() == -0)
            p.N.z() = 0.0;

        // B_i = (T_i * N_i).normalized
        p.B = (Vector3f::cross(p.T, p.N)).normalized();

        // get rid of -0 values
        if(p.B.z() == -0)
            p.B.z() = 0.0;

        // store binormal
        B_prev = p.B;

        /*cout << "V = (" <<   p.V.x() << ", " <<    p.V.y() << ", " <<
        p.V.z() << ")" << endl;

        cout << "T = (" <<   p.T.x() << ", " <<    p.T.y() << ", " <<
        p.T.z() << ")" << endl;

        cout << "N = (" <<   p.N.x() << ", " <<    p.N.y() << ", " <<
        p.N.z() << ")" << endl;

        cout << "B = (" <<   p.B.x() << ", " <<    p.B.y() << ", " <<
        p.B.z() << ")" << endl;*/

        // add point to curve
        curve.push_back(p);
    }

    return curve;
}








Curve evalBspline( const vector< Vector3f >& P, unsigned steps )
{
    // Check
    if( P.size() < 4 )
    {
        cerr << "evalBspline must be called with 4 or more control points." << endl;
        exit( 0 );
    }

    // curve we will be returning
    Curve finalCurve;

    // segment of final curve
    Curve curveSegment;

    // loop through control points after the initial 4
    for(unsigned i = 3; i < P.size(); i++)
    {
        // store the new control point i in last col of geometric matrix
        Matrix4f geo_matrix(
                            P[i-3].x(), P[i-2].x(), P[i-1].x(), P[i].x(),
                            P[i-3].y(), P[i-2].y(), P[i-1].y(), P[i].y(),
                            P[i-3].z(), P[i-2].z(), P[i-1].z(), P[i].z(),
                            0.0f, 0.0f, 0.0f, 1.0f
                            );

        //cout << "current control points\n";

        //geo_matrix.print();

        // vector to store our transformed control points
        vector< Vector3f > new_cPoints;

        // get new control points by changing from b_spline to bernstein basis
        // we can then pass them to the evalBezier function
        Matrix4f new_geo_matrix = geo_matrix * b_spline * bernstein.inverse();

        //cout << "\nnew geo matrix\n";

        //new_geo_matrix.print();

        // for each of the 4 new control points
        for(int k = 0; k < 4; k++)
        {
            // get cp coords from the columns of the matrix
            Vector3f c_point(new_geo_matrix(0,k), new_geo_matrix(1,k), new_geo_matrix(2,k));

            //c_point.print();

            new_cPoints.push_back(c_point);
        }

        //cout << "\nnew control points\n";

        //for(int d = 0; d < 4; d++)
            //cout << "(" << new_cPoints[d].x() << ", " << new_cPoints[d].y() << ", " << new_cPoints[d].z() << ")\n";

        // pass control points to evalBezier to get our curve segment
        Curve curveSegment = evalBezier(new_cPoints, steps);

        // add curve segment to our curve
        for(unsigned j = 0; j < curveSegment.size(); j++)
        {
            CurvePoint cp = curveSegment[j];

            finalCurve.push_back(cp);
        }
    }

    return finalCurve;
}








Curve evalCircle( float radius, unsigned steps )
{
    // This is a sample function on how to properly initialize a Curve
    // (which is a vector< CurvePoint >).
    
    // Preallocate a curve with steps+1 CurvePoints
    Curve R( steps+1 );

    // Fill it in counterclockwise
    for( unsigned i = 0; i <= steps; ++i )
    {
        // step from 0 to 2pi
        float t = 2.0f * M_PI * float( i ) / steps;

        // Initialize position
        // We're pivoting counterclockwise around the y-axis
        R[i].V = radius * Vector3f( cos(t), sin(t), 0 );
        
        // Tangent vector is first derivative
        R[i].T = Vector3f( -sin(t), cos(t), 0 );
        
        // Normal vector is second derivative
        R[i].N = Vector3f( -cos(t), -sin(t), 0 );

        // Finally, binormal is facing up.
        R[i].B = Vector3f( 0, 0, 1 );
    }

    return R;
}

void drawCurve( const Curve& curve, float framesize )
{
    // Save current state of OpenGL
    glPushAttrib( GL_ALL_ATTRIB_BITS );

    // Setup for line drawing
    glDisable( GL_LIGHTING ); 
    glColor4f( 1, 1, 1, 1 );
    glLineWidth( 1 );
    
    // Draw curve
    glBegin( GL_LINE_STRIP );
    for( unsigned i = 0; i < curve.size(); ++i )
    {
        glVertex( curve[ i ].V );
    }
    glEnd();

    glLineWidth( 1 );

    // Draw coordinate frames if framesize nonzero
    if( framesize != 0.0f )
    {
        Matrix4f M;

        for( unsigned i = 0; i < curve.size(); ++i )
        {
            M.setCol( 0, Vector4f( curve[i].N, 0 ) );
            M.setCol( 1, Vector4f( curve[i].B, 0 ) );
            M.setCol( 2, Vector4f( curve[i].T, 0 ) );
            M.setCol( 3, Vector4f( curve[i].V, 1 ) );

            glPushMatrix();
            glMultMatrixf( M );
            glScaled( framesize, framesize, framesize );
            glBegin( GL_LINES );
            glColor3f( 1, 0, 0 ); glVertex3d( 0, 0, 0 ); glVertex3d( 1, 0, 0 );
            glColor3f( 0, 1, 0 ); glVertex3d( 0, 0, 0 ); glVertex3d( 0, 1, 0 );
            glColor3f( 0, 0, 1 ); glVertex3d( 0, 0, 0 ); glVertex3d( 0, 0, 1 );
            glEnd();
            glPopMatrix();
        }
    }
    
    // Pop state
    glPopAttrib();
}

