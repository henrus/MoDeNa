/**
@cond

   ooo        ooooo           oooooooooo.             ooooo      ooo
   `88.       .888'           `888'   `Y8b            `888b.     `8'
    888b     d'888   .ooooo.   888      888  .ooooo.   8 `88b.    8   .oooo.
    8 Y88. .P  888  d88' `88b  888      888 d88' `88b  8   `88b.  8  `P  )88b
    8  `888'   888  888   888  888      888 888ooo888  8     `88b.8   .oP"888
    8    Y     888  888   888  888     d88' 888    .o  8       `888  d8(  888
   o8o        o888o `Y8bod8P' o888bood8P'   `Y8bod8P' o8o        `8  `Y888""8o

Copyright
    2014-2016 MoDeNa Consortium, All rights reserved.

License
    This file is part of Modena.

    The Modena interface library is free software; you can redistribute it
    and/or modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    Modena is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along
    with Modena.  If not, see <http://www.gnu.org/licenses/>.

@endcond
@file
Low-level interface library
@author    Henrik Rusche
@copyright 2014-2016, MoDeNa Project. GNU Public License.
@defgroup  C++_interface_library
MoDeNa C++ interface library
*/

#ifndef modenScalarField_H
#define modenScalarField_H

#include "volFields.H"
#include "DynamicList.H"
#include "modenaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class errorOp
{
public:

    int operator()(const int& x, const int& y) const
    {
        if(y > 0)
        {
            return y;
        }
        
        return x;
    }
};


class modenaScalarField
{
    // Private data

        // MoDeNa model
        Modena::modenaModel model_;

        class iTuple
        {
            size_t pos_;
            const volScalarField* field_;

        public:
            iTuple()
            :
                pos_(0),
                field_(NULL)
            {
            }

            iTuple(size_t pos, const volScalarField& field)
            :
                pos_(pos),
                field_(&field)
            {
            }

            size_t pos() const
            {
                return pos_;
            }

            const volScalarField& field() const
            {
                return *field_;
            }
        };

        class oTuple
        {
            size_t pos_;
            volScalarField* field_;

        public:
            oTuple()
            :
                pos_(0),
                field_(NULL)
            {
            }

            oTuple(size_t pos, volScalarField& field)
            :
                pos_(pos),
                field_(&field)
            {
            }

            size_t pos() const
            {
                return pos_;
            }
            
            volScalarField& field()
            {
                return *field_;
            }
        };

        DynamicList<iTuple> inputs_;
        DynamicList<oTuple> outputs_;

        void applyModel()
        {
            const volScalarField& l = inputs_[0].field();
            forAll(l, cellI)
            {
                forAll(inputs_, fieldI)
                {
                    model_.inputs_set
                    (
                        inputs_[fieldI].pos(),
                        inputs_[fieldI].field()[cellI]
                    );
                }

                model_.call();

                forAll(outputs_, fieldI)
                {
                    outputs_[fieldI].field()[cellI] =
                        model_.outputs_get(outputs_[fieldI].pos());
                }
            }

            forAll(l.boundaryField(), boundaryI)
            {
                forAll(l.boundaryField()[boundaryI], faceI)
                {
                    forAll(inputs_, fieldI)
                    {
                        const iTuple& t = inputs_[fieldI];
                        const fvPatchScalarField& pf =
                            t.field().boundaryField()[boundaryI];
                        model_.inputs_set(t.pos(), pf[faceI]);
                    }

                    model_.call();

                    forAll(outputs_, fieldI)
                    {
                        oTuple& t = outputs_[fieldI];
                        t.field().boundaryFieldRef()[boundaryI][faceI] =
                            model_.outputs_get(t.pos());
                    }
                }
            }

            forAll(outputs_, fieldI)
            {
                outputs_[fieldI].field().correctBoundaryConditions();
            }
        }


public:
 
    // Constructors

        //- Construct from model name
        modenaScalarField(const word& name)
        try:
            model_(name),
            inputs_(),
            outputs_()
        {
        }
        catch(const Modena::modenaException& e)
        {
            std::exit(e.errorCode());
        }

    // Destructor

        ~modenaScalarField(){};


    // Member Functions

        Modena::modenaModel& model()
        {
            return model_;
        }
    
        const Modena::modenaModel& model() const
        {
            return model_;
        }

        void registerInputField(const word& name, const volScalarField& field)
        {
            inputs_.append(iTuple(model_.inputs_argPos(name), field));
        }

        void registerOutputField(const word& name, volScalarField& field)
        {
            outputs_.append(oTuple(model_.inputs_argPos(name), field));
        }

        void argPosCheck() const
        {
            model_.argPos_check();
        }

        void callAndExit()
        {
            try
            {
                applyModel();
            }
            catch(const Modena::modenaException& e)
            {
                std::exit(e.errorCode());
            }
        }

        void call()
        {
            try
            {
                applyModel();
            }
            catch(const Modena::modenaException& e)
            {
                throw Modena::modenaException(e.errorCode());
            }
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
