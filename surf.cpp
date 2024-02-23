#include "surf.h"
#include "extra.h"
using namespace std;

namespace
{
    
    // We're only implenting swept surfaces where the profile curve is
    // flat on the xy-plane.  This is a check function.
    static bool checkFlat(const Curve &profile)
    {
        for (unsigned i=0; i<profile.size(); i++)
            if (profile[i].V[2] != 0.0 ||
                profile[i].T[2] != 0.0 ||
                profile[i].N[2] != 0.0)
                return false;
    
        return true;
    }
}






Surface makeSurfRev(const Curve &profile, unsigned steps)
{
    Surface surface;
    
    if (!checkFlat(profile))
    {
        cerr << "surfRev profile curve must be flat on xy plane." << endl;
        exit(0);
    }

    // Implementing the surfaces of revolution
    // Step 1: generate the vertices of the points of revolution
    for(size_t i = 0; i < profile.size(); i++)
    {

        // Now, rotate each point of the curve around the y axis, then add to vertices
        // This involves calculating the new location of the point based on the angle of change
        for(unsigned j = 0; j <= steps; j++)  
        {
            // Angle that determines how many degrees we will rotate
            float angle = 2 * M_PI * j / steps;

            // Get the current point we are working with
            Vector4f profilePoint = Vector4f(profile[i].V[0], profile[i].V[1], 
            profile[i].V[2], 1.0f); 

            // Your matrix setup is fine
            Matrix4f rotateM = rotateM.rotateY(angle);
            Matrix3f subrotateM = rotateM.getSubmatrix3x3(0,0);
            Matrix3f tsubrotateM = subrotateM.transposed();
            Matrix3f rotateNormalM = tsubrotateM.inverse();

            //Multiply the profile point with the rotate y matrix and save in a Vector3f
            Vector4f profilePointRotated= rotateM * profilePoint ;
            Vector3f surfaceVV = Vector3f(profilePointRotated[0], profilePointRotated[1], 
            profilePointRotated[2]);

            // Multiply the normal matrix with the profilePoint normal.
            Vector3f surfaceVN = rotateNormalM * profile[i].N; 

            // Push to Curve list
            surface.VV.push_back(surfaceVV );
            surface.VN.push_back(-1 * surfaceVN);
            
        }
    }

    // generate faces
    for(unsigned i = 0; i < surface.VV.size() - steps + 1; i++) {
        // To hold the current two triangles that make up the current square
        Tup3u triangle1;
        Tup3u triangle2;

        if((i + 1) % steps + 1 != 0) {
            triangle1 = Tup3u(i + 1, i, i + steps + 1);
            triangle2 = Tup3u(i, i + steps, i + steps + 1);
        }

        if(i == 1 || i == 2) {
            cout << triangle1 << endl;
            cout << triangle2 << endl;
        }
        surface.VF.push_back(triangle1);
        surface.VF.push_back(triangle2);
        
    }



    
    cout << "right before returning the surface" << endl;
    return surface;
}




Surface makeGenCyl(const Curve &profile, const Curve &sweep )
{
    Surface surface;

    if (!checkFlat(profile))
    {
        cerr << "genCyl profile curve must be flat on xy plane." << endl;
        exit(0);
    }

    // TODO: Here you should build the surface.  See surf.h for details.

    for (unsigned i = 0; i < profile.size(); i++) {
        for (unsigned j = 0; j < sweep.size(); j++) {
            // Matrix for coordinate transformation
            Matrix4f coords_M
            (   sweep[j].N[0], sweep[j].B[0], sweep[j].T[0], sweep[j].V[0],
                sweep[j].N[1], sweep[j].B[1], sweep[j].T[1], sweep[j].V[1],
                sweep[j].N[2], sweep[j].B[2], sweep[j].T[2], sweep[j].V[2],
                    0.0f,          0.0f,          0.0f,           1.0f     );

            // unsure how to calculate angle for this part
            float angle = 2 * M_PI * j / profile.size();
            // Calculate the matrix we need for determining the new normals
            Matrix4f rotateM = Matrix4f::rotateY(angle);

            Matrix3f subrotateM = rotateM.getSubmatrix3x3(0,0);
            Matrix3f tSubrotateM = subrotateM.transposed();
            Matrix3f rotateNormalM = tSubrotateM.inverse();

            // Calculating the new vertex of the coordinates around the y axis
            Vector4f profileVV = Vector4f(profile[i].V[0], profile[i].V[1], profile[i].V[2], 1.0f);
            Vector4f coordsVV = coords_M * profileVV;
            Vector3f surfaceVV = Vector3f(coordsVV[0], coordsVV[1], coordsVV[2]);

            // Now that we have the new vertex, we need to calculate its normal as well
            Vector3f surfaceVN = rotateNormalM * profile[i].N * -1;

            // Now that we have the new vertices and their normals, add them to the surface object
            surface.VV.push_back(surfaceVV);
            surface.VN.push_back(surfaceVN);

        }
    }

    // Now, generate the faces
    for (unsigned k = 0; k < profile.size() - 1; k++) {
        for (unsigned j = 0; j < sweep.size() - 1; j++) {
        // Define the indices of the vertices for the current face
        unsigned v0 = k * sweep.size() + j;
        unsigned v1 = v0 + 1;
        unsigned v2 = (k + 1) * sweep.size() + j;
        unsigned v3 = v2 + 1;

        // Add the indices to form two triangles
        surface.VF.push_back(Tup3u(v0, v2, v1));
        surface.VF.push_back(Tup3u(v1, v2, v3));
    }
}



    cerr << "\t>>> makeGenCyl called (but not implemented).\n\t>>> Returning empty surface." <<endl;

    return surface;
}





void drawSurface(const Surface &surface, bool shaded)
{
    // Save current state of OpenGL
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    if (shaded)
    {
        // This will use the current material color and light
        // positions.  Just set these in drawScene();
        glEnable(GL_LIGHTING);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        // This tells openGL to *not* draw backwards-facing triangles.
        // This is more efficient, and in addition it will help you
        // make sure that your triangles are drawn in the right order.
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
    }
    else
    {        
        glDisable(GL_LIGHTING);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        
        glColor4f(0.4f,0.4f,0.4f,1.f);
        glLineWidth(1);
    }

    glBegin(GL_TRIANGLES);
    for (unsigned i=0; i<surface.VF.size(); i++)
    {
        glNormal(surface.VN[surface.VF[i][0]]);
        glVertex(surface.VV[surface.VF[i][0]]);
        glNormal(surface.VN[surface.VF[i][1]]);
        glVertex(surface.VV[surface.VF[i][1]]);
        glNormal(surface.VN[surface.VF[i][2]]);
        glVertex(surface.VV[surface.VF[i][2]]);
    }
    glEnd();

    glPopAttrib();
}

void drawNormals(const Surface &surface, float len)
{
    // Save current state of OpenGL
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    glDisable(GL_LIGHTING);
    glColor4f(0,1,1,1);
    glLineWidth(1);

    glBegin(GL_LINES);
    for (unsigned i=0; i<surface.VV.size(); i++)
    {
        glVertex(surface.VV[i]);
        glVertex(surface.VV[i] + surface.VN[i] * len);
    }
    glEnd();

    glPopAttrib();
}

void outputObjFile(ostream &out, const Surface &surface)
{
    
    for (unsigned i=0; i<surface.VV.size(); i++)
        out << "v  "
            << surface.VV[i][0] << " "
            << surface.VV[i][1] << " "
            << surface.VV[i][2] << endl;

    for (unsigned i=0; i<surface.VN.size(); i++)
        out << "vn "
            << surface.VN[i][0] << " "
            << surface.VN[i][1] << " "
            << surface.VN[i][2] << endl;

    out << "vt  0 0 0" << endl;
    
    for (unsigned i=0; i<surface.VF.size(); i++)
    {
        out << "f  ";
        for (unsigned j=0; j<3; j++)
        {
            unsigned a = surface.VF[i][j]+1;
            out << a << "/" << "1" << "/" << a << " ";
        }
        out << endl;
    }
}
