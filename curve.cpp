#include "curve.h"
#include "extra.h"
#ifdef WIN32
#include <windows.h>
#endif
#include <GL/gl.h>
#include <fstream>
using namespace std;

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

    ofstream out;

    out.open("curve_points.dat");

    Curve curve;

    float t;

    Vector3f B_prev(1, 1, 1);

    Matrix4f bernstein(
                            1.0f, -3.0f, 3.0f, -1.0f,
                            0.0f, 3.0f, -6.0f, 3.0f,
                            0.0f, 0.0f, 3.0f, -3.0f,
                            0.0f, 0.0f, 0.0f, 1.0f
                        );


    Matrix4f d_bernstein(
                            0.0f, -3.0f, 6.0f, -3.0f,
                            0.0f, 3.0f, -12.0f, 9.0f,
                            0.0f, 0.0f, 6.0f, -9.0f,
                            0.0f, 0.0f, 0.0f, 3.0f
                            );

    Matrix4f geo_matrix(
                            P[0].x(), P[1].x(), P[2].x(), P[3].x(),
                            P[0].y(), P[1].y(), P[2].y(), P[3].y(),
                            P[0].z(), P[1].z(), P[2].z(), P[3].z(),
                            0.0f, 0.0f, 0.0f, 1.0f
                            );

    

    // Determine the number of points to generate on this segment based on the steps parameter
    int num_points = steps + 1;

    out << "\t>>> evalBezier has been called with the following input:" << endl;

    out << "\t>>> Control points (type vector< Vector3f >): "<< endl;
    for( unsigned i = 0; i < P.size(); ++i )
    {
        out << "\t>>> (" << P[i].x() << ", " << P[i].y() << ", " << P[i].z() << ")" << endl;
    }

    out << "\t>>> Steps (type steps): " << steps << endl;


    for(int k = 0; k < num_points - 1; k++) 
    {
        CurvePoint p;

        Vector4f d_q_t;

        Vector4f T;
        
        t = static_cast<float>(k) / static_cast<float>(steps);

        Vector4f monomials(1, t, pow(t, 2), pow(t, 3));
        
        out << "the value of t is: " << t << endl;

        Vector4f result = geo_matrix * bernstein * monomials;

        p.V = Vector3f(result.x(), result.y(), result.z());

        d_q_t = geo_matrix * d_bernstein * monomials;

        T = d_q_t.normalized();

        p.T = Vector3f(T.x(), T.y(), T.z());

        p.N = (B_prev * p.T).normalized();

        p.B = (p.N * p.T).normalized();

        out << "(" <<   p.V.x() << ", " <<    p.V.y() << ", " <<
        p.V.z() << ")" << endl;

        out << "(" <<   p.T.x() << ", " <<    p.T.y() << ", " <<
        p.T.z() << ")" << endl;

        out << "(" <<   p.N.x() << ", " <<    p.N.y() << ", " <<
        p.N.z() << ")" << endl;

        out << "(" <<   p.B.x() << ", " <<    p.B.y() << ", " <<
        p.B.z() << ")" << endl;

        curve.push_back(p);

        B_prev = p.B;
    }

    out.close();

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

    // TODO:
    // It is suggested that you implement this function by changing
    // basis from B-spline to Bezier.  That way, you can just call
    // your evalBezier function.

    // new_geo_matrix = geo_matrix * b_spline * inverse_bernstein

    // evalBezier(new_geo_matrix, steps)

    cerr << "\t>>> evalBSpline has been called with the following input:" << endl;

    cerr << "\t>>> Control points (type vector< Vector3f >): "<< endl;
    for( unsigned i = 0; i < P.size(); ++i )
    {
        cerr << "\t>>> " << P[i] << endl;
    }

    cerr << "\t>>> Steps (type steps): " << steps << endl;
    cerr << "\t>>> Returning empty curve." << endl;

    // Return an empty curve right now.
    return Curve();
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

