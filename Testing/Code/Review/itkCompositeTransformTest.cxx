/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkCompositeTransformTest.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <iostream>

#include "itkAffineTransform.h"
#include "itkScaleTransform.h"
#include "itkCompositeTransform.h"
#include "itkArray2D.h"

namespace {

const double epsilon = 1e-10;

template <typename TPoint>
bool testPoint( const TPoint & p1, const TPoint & p2 )
  {
  bool pass=true;
  for ( unsigned int i = 0; i < TPoint::PointDimension; i++ )
    {
    if( vcl_fabs( p1[i] - p2[i] ) > epsilon )
      pass=false;
    }
  return pass;
  }

template <typename TMatrix>
bool testMatrix( const TMatrix & m1, const TMatrix & m2 )
  {
  unsigned int i, j;
  bool pass=true;
  for ( i = 0; i < TMatrix::RowDimensions; i++ )
    {
    for ( j = 0; j < TMatrix::ColumnDimensions; j++ )
      {
      if( vcl_fabs( m1[i][j] - m2[i][j] ) > epsilon )
        pass=false;
      }
  }
  return pass;
  }

template <typename TArray2D>
bool testJacobian( const TArray2D & m1, const TArray2D & m2 )
  {
  unsigned int i, j;
  bool pass=true;
  for ( i = 0; i < m1.rows(); i++ )
    {
    for ( j = 0; j < m1.cols(); j++ )
      {
      if( vcl_fabs( m1[i][j] - m2[i][j] ) > epsilon )
        pass=false;
      }
  }
  return pass;
  }

template <typename TVector>
bool testVectorArray( const TVector & v1, const TVector & v2 )
  {
  bool pass=true;
  for ( unsigned int i = 0; i < v1.Size(); i++ )
    {
    if( vcl_fabs( v1[i] - v2[i] ) > epsilon )
      pass=false;
    }
  return pass;
  }

} //namespace

/******/

int itkCompositeTransformTest(int ,char *[] )
{
  const unsigned int NDimensions = 2;

  /* Create composite transform */
  typedef itk::CompositeTransform<double, NDimensions>  CompositeType;
  typedef CompositeType::ScalarType                     ScalarType;

  CompositeType::Pointer compositeTransform = CompositeType::New();

  /* Test obects */
  typedef  itk::Matrix<ScalarType,NDimensions,NDimensions>   Matrix2Type;
  typedef  itk::Vector<ScalarType,NDimensions>               Vector2Type;

  /* Test that we have an empty the queue */
  if( compositeTransform->GetNumberOfSubTransforms() != 0 )
    {
    std::cout << "Failed. Expected GetNumberOfSubTransforms to return 0."
              << std::endl;
    return EXIT_FAILURE;
    }

  /* Add an affine transform */
  typedef itk::AffineTransform<ScalarType, NDimensions> AffineType;
  AffineType::Pointer affine = AffineType::New();
  Matrix2Type matrix2;
  Vector2Type vector2;
  matrix2[0][0] = 1;
  matrix2[0][1] = 2;
  matrix2[1][0] = 3;
  matrix2[1][1] = 4;
  vector2[0] = 5;
  vector2[1] = 6;
  affine->SetMatrix(matrix2);
  affine->SetOffset(vector2);

  compositeTransform->AddTransform( affine );

  /* Test that we have one tranform in the queue */
  if( compositeTransform->GetNumberOfSubTransforms() != 1 )
    {
    std::cout << "Failed adding transform to queue." << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Composite Transform:" << std::endl << compositeTransform;

  /* Retrieve the transform and check that it's the same */
  AffineType::ConstPointer affineGet;
  affineGet = dynamic_cast<AffineType const *>
    ( compositeTransform->GetNthTransform(0).GetPointer() );
  if( affineGet.IsNull() )
    {
    std::cout << "Failed retrieving transform from queue." << std::endl;
    return EXIT_FAILURE;
    }

  Matrix2Type matrix2Get = affineGet->GetMatrix();
  Vector2Type vector2Get = affineGet->GetOffset();
  if( !testMatrix(matrix2, matrix2Get) || !testVectorArray(vector2, vector2Get ) )
    {
    std::cout << "Failed retrieving correct transform data." << std::endl;
    return EXIT_FAILURE;
    }

  /* Get parameters with single transform.
   * Should be same as GetParameters from affine transform. */
  CompositeType::ParametersType parametersTest, parametersTruth;
  parametersTest = compositeTransform->GetParameters();
  parametersTruth = affine->GetParameters();
  std::cout << "Get Parameters: " << std::endl
            << "affine parametersTruth: " << std::endl << parametersTruth
            << std::endl
            << "parametersTest from Composite: " << std::endl << parametersTest
            << std::endl;

  if( !testVectorArray( parametersTest, parametersTruth ) )
    {
    std::cout << "Failed GetParameters() for single transform." << std::endl;
    return EXIT_FAILURE;
    }

  /* Set parameters with single transform. */
  CompositeType::ParametersType parametersNew(6), parametersReturned;
  parametersNew[0] = 0;
  parametersNew[1] = 10;
  parametersNew[2] = 20;
  parametersNew[3] = 30;
  parametersNew[4] = 40;
  parametersNew[5] = 50;
  std::cout << "Set Parameters: " << std::endl;
  compositeTransform->SetParameters( parametersNew );
  std::cout << "retrieving... " << std::endl;
  parametersReturned = compositeTransform->GetParameters();
  std::cout << "parametersNew: " << std::endl << parametersNew << std::endl
            << "parametersReturned: " << std::endl << parametersReturned
            << std::endl;
  if( !testVectorArray( parametersNew, parametersReturned ) )
    {
    std::cout << "Failed SetParameters() for single transform." << std::endl;
    return EXIT_FAILURE;
    }

  /* Test fixed parameters set/get */
  parametersTest = compositeTransform->GetFixedParameters();
  parametersTruth = affine->GetFixedParameters();
  std::cout << "Get Fixed Parameters: " << std::endl
            << "affine parametersTruth: " << std::endl << parametersTruth
            << std::endl
            << "parametersTest from Composite: " << std::endl << parametersTest
            << std::endl;

  if( !testVectorArray( parametersTest, parametersTruth ) )
    {
    std::cout << "Failed GetFixedParameters() for single transform." << std::endl;
    return EXIT_FAILURE;
    }

  parametersNew.SetSize( parametersTruth.Size() );
  parametersNew.Fill(1);
  parametersNew[0] = 42;

  std::cout << "Set Fixed Parameters: " << std::endl;
  compositeTransform->SetFixedParameters( parametersNew );
  std::cout << "retrieving... " << std::endl;
  parametersReturned = compositeTransform->GetFixedParameters();
  std::cout << "parametersNew: " << std::endl << parametersNew << std::endl
            << "parametersReturned: " << std::endl << parametersReturned
            << std::endl;
  if( !testVectorArray( parametersNew, parametersReturned ) )
    {
    std::cout << "Failed SetFixedParameters() for single transform." << std::endl;
    return EXIT_FAILURE;
    }

  /* Reset affine transform to original values */
  compositeTransform->ClearTransformQueue();
  affine = AffineType::New();
  affine->SetMatrix(matrix2);
  affine->SetOffset(vector2);
  compositeTransform->AddTransform( affine );

  /* Setup test point and truth value for tests */
  CompositeType::InputPointType  inputPoint;
  CompositeType::OutputPointType outputPoint, affineTruth;
  inputPoint[0] = 2;
  inputPoint[1] = 3;
  affineTruth[0] = 13;
  affineTruth[1] = 24;

  /* Test transforming the point with just the single affine transform */
  outputPoint = compositeTransform->TransformPoint( inputPoint );
  if( !testPoint( outputPoint, affineTruth) )
      {
      std::cout << "Failed transforming point with single transform."
                << std::endl;
      return EXIT_FAILURE;
      }

  /* Test inverse */
  CompositeType::Pointer inverseTransform = CompositeType::New();
  CompositeType::OutputPointType inverseTruth, inverseOutput;
  if( !compositeTransform->GetInverse( inverseTransform ) )
    {
    std::cout << "Expected GetInverse() to succeed." << std::endl;
    return EXIT_FAILURE;
    }
  inverseTruth = inputPoint;
  inverseOutput = inverseTransform->TransformPoint( affineTruth );
  std::cout << "Transform point with inverse composite transform: "
            << std::endl
            << "  Test point: " << affineTruth << std::endl
            << "  Truth: " << inverseTruth << std::endl
            << "  Output: " << inverseOutput << std::endl;
  if( !testPoint( inverseOutput, inverseTruth ) )
    {
    std::cout << "Failed transform point with inverse composite transform (1)."
              << std::endl;
    return EXIT_FAILURE;
    }

  /* Test GetJacobian */
  CompositeType::JacobianType jacComposite, jacSingle;
  CompositeType::InputPointType jacPoint;
  jacPoint[0]=1;
  jacPoint[1]=2;
  jacComposite = compositeTransform->GetJacobian( jacPoint );
  jacSingle = affine->GetJacobian( jacPoint );
  std::cout << "Single jacobian:" << std::endl << jacSingle << std::endl;
  if( !testJacobian( jacComposite, jacSingle ) )
    {
    std::cout << "Failed getting jacobian for single transform." << std::endl;
    return EXIT_FAILURE;
    }

  /*
   * Create and add 2nd transform
   */

  typedef itk::ScaleTransform<double, NDimensions>  ScaleTransformType;
  ScaleTransformType::Pointer scaleTransform = ScaleTransformType::New();

  ScaleTransformType::ScaleType::ValueType iscaleInit[2] = {2,-0.5};
  ScaleTransformType::ScaleType            iscale = iscaleInit;
  scaleTransform->SetScale( iscale );

  compositeTransform->AddTransform( scaleTransform );

  std::cout << std::endl << "Two-component Composite Transform:"
            << std::endl << compositeTransform;
  std::cout << std::endl << "Transform at position 0: "
            << std::endl << compositeTransform->GetNthTransform( 0 );

  /* Test that we have two tranforms in the queue */
  if( compositeTransform->GetNumberOfSubTransforms() != 2 )
    {
    std::cout << "Failed adding 2nd transform to queue." << std::endl;
    return EXIT_FAILURE;
    }

  /* Transform a point with both transforms */
  CompositeType::OutputPointType compositeTruth;
  compositeTruth[0] = affineTruth[0] * iscale[0];
  compositeTruth[1] = affineTruth[1] * iscale[1];

  outputPoint = compositeTransform->TransformPoint( inputPoint );
  std::cout << "Transform point with two-component composite transform: "
            << std::endl
            << "  Test point: " << inputPoint << std::endl
            << "  Truth: " << compositeTruth << std::endl
            << "  Output: " << outputPoint << std::endl;

  if( !testPoint( outputPoint, compositeTruth) )
    {
    std::cout << "Failed transforming point with two transforms."
              << std::endl;
    return EXIT_FAILURE;
    }

  /* Test inverse with two transforms */
  if( !compositeTransform->GetInverse( inverseTransform ) )
    {
    std::cout << "Expected GetInverse() to succeed." << std::endl;
    return EXIT_FAILURE;
    }
  inverseTruth = inputPoint;
  inverseOutput = inverseTransform->TransformPoint( compositeTruth );
  std::cout << "Transform point with two-component inverse composite transform: "
            << std::endl
            << "  Test point: " << compositeTruth << std::endl
            << "  Truth: " << inverseTruth << std::endl
            << "  Output: " << inverseOutput << std::endl;
  if( !testPoint( inverseOutput, inverseTruth ) )
    {
    std::cout << "Failed transform point with two-component inverse composite transform."
              << std::endl;
    return EXIT_FAILURE;
    }

  /* Get inverse transform again, but using other accessor. */
  CompositeType::ConstPointer inverseTransform2;
  inverseTransform2 = dynamic_cast<const CompositeType *>
    ( compositeTransform->GetInverseTransform().GetPointer() );
  if( ! inverseTransform2 )
    {
    std::cout << "Failed calling GetInverseTransform()." << std::endl;
    return EXIT_FAILURE;
    }
  inverseOutput = inverseTransform2->TransformPoint( compositeTruth );
  if( !testPoint( inverseOutput, inverseTruth ) )
    {
    std::cout << "Failed transform point with two-component inverse composite transform (2)."
              << std::endl;
    return EXIT_FAILURE;
    }

  /* Test IsLinear() by calling on each sub transform */
  bool allAreLinear = true;
  for( unsigned int n=0;
        n < compositeTransform->GetNumberOfSubTransforms(); n ++)
    {
    if( ! compositeTransform->GetNthTransform( n )->IsLinear() )
      {
      allAreLinear = false;
      }
    }
  if( compositeTransform->IsLinear() != allAreLinear )
    {
    std::cout << "compositeTransform returned unexpected value for IsLinear()."
      " Expected " << allAreLinear << std::endl;
    return EXIT_FAILURE;
    }

  /* Get parameters with multi-transform. */
  parametersTest = compositeTransform->GetParameters();
  unsigned int affineParamsN = affine->GetNumberOfParameters();
  unsigned int scaleParamsN = scaleTransform->GetNumberOfParameters();
  parametersTruth.SetSize( scaleParamsN + affineParamsN );
  /* Fill using different method than is used in the class.
     Remember we added scaleTransform 2nd, so it's at front of queue */
  for( unsigned int n=0; n < scaleParamsN; n++)
    parametersTruth.SetElement(
      n, scaleTransform->GetParameters().GetElement( n ) );
  for( unsigned int n=0; n < affineParamsN; n++)
    parametersTruth.SetElement( n + scaleParamsN,
      affine->GetParameters().GetElement( n ) );
  std::cout << "Get Multi-transform Parameters: " << std::endl
            << "parametersTruth: " << std::endl << parametersTruth
            << std::endl
            << "parametersTest from Composite: " << std::endl << parametersTest
            << std::endl;

  if( !testVectorArray( parametersTest, parametersTruth ) )
    {
    std::cout << "Failed GetParameters() for multi transform." << std::endl;
    return EXIT_FAILURE;
    }

  /* Set parameters with multi transform. */
  parametersNew.SetSize( parametersTruth.Size() );
  parametersNew.Fill( -1 );
  parametersNew[0] = 19;
  parametersNew[ parametersTruth.Size() - 1 ] = 71;
  std::cout << "Set Multi-transform Parameters: " << std::endl;
  compositeTransform->SetParameters( parametersNew );
  std::cout << "retrieving... " << std::endl;
  parametersReturned = compositeTransform->GetParameters();
  std::cout << "parametersNew: " << std::endl << parametersNew << std::endl
            << "parametersReturned: " << std::endl << parametersReturned
            << std::endl;
  if( !testVectorArray( parametersNew, parametersReturned ) )
    {
    std::cout << "Failed SetParameters() for multi transform." << std::endl;
    return EXIT_FAILURE;
    }

  /* Test get fixed parameters with multi-transform */
  parametersTest = compositeTransform->GetFixedParameters();
  affineParamsN = affine->GetFixedParameters().Size();
  scaleParamsN = scaleTransform->GetFixedParameters().Size();
  parametersTruth.SetSize( scaleParamsN + affineParamsN );
  for( unsigned int n=0; n < scaleParamsN; n++)
    parametersTruth.SetElement(
      n, scaleTransform->GetFixedParameters().GetElement( n ) );
  for( unsigned int n=0; n < affineParamsN; n++)
    parametersTruth.SetElement( n + scaleParamsN,
      affine->GetFixedParameters().GetElement( n ) );
  std::cout << "Get Multi-transform Fixed Parameters: " << std::endl
            << "parametersTruth: " << std::endl << parametersTruth
            << std::endl
            << "parametersTest: " << std::endl << parametersTest
            << std::endl;

  if( !testVectorArray( parametersTest, parametersTruth ) )
    {
    std::cout << "Failed GetFixedParameters() for multi transform." << std::endl;
    return EXIT_FAILURE;
    }

  /* Test set fixed parameters with multi-transform */
  /* NOTE that ScaleTransform doesn't have any fixed transforms. 'Get' returns
   * something, but 'Set' does nothing, just an empty function */

  parametersNew = parametersTruth;
  /* Leave the ScaleTransform part of this as-is since the set func
   * doesn't do anything. */
  parametersNew[ scaleParamsN ] = 123;
  parametersNew[ parametersTruth.Size() - 1 ] = 246;

  std::cout << "Set Multi-transform Fixed Parameters: " << std::endl;
  compositeTransform->SetFixedParameters( parametersNew );
  std::cout << "retrieving... " << std::endl;
  parametersReturned = compositeTransform->GetFixedParameters();
  std::cout << "parametersNew: " << std::endl << parametersNew << std::endl
            << "parametersReturned: " << std::endl << parametersReturned
            << std::endl;
  //std::cout << "Composite Transform: " << std::endl << compositeTransform;
  if( !testVectorArray( parametersNew, parametersReturned ) )
    {
    std::cout << "Failed SetFixedParameters() for multi transform." << std::endl;
    return EXIT_FAILURE;
    }

  unsigned int nParameters = compositeTransform->GetNumberOfParameters();
  std::cout << "Number of parameters: " << nParameters << std::endl;

  /* Test GetJacobian with two transforms */
  CompositeType::JacobianType jacTruth, jacComposite2, jacAffine, jacScale;
  CompositeType::InputPointType jacPoint2;
  jacPoint2[0]=1;
  jacPoint2[1]=2;
  jacComposite2 = compositeTransform->GetJacobian( jacPoint2 );
  jacAffine = affine->GetJacobian( jacPoint2 );
  jacScale = scaleTransform->GetJacobian( jacPoint2 );
  jacTruth.SetSize( jacComposite2.rows(), jacComposite2.cols() );
  jacTruth.update( jacScale, 0, 0 );
  jacTruth.update( jacAffine, 0, jacScale.cols() );
  std::cout << "Affine jacobian:" << std::endl << jacAffine;
  std::cout << "Scale jacobian:" << std::endl << jacScale;
  std::cout << "Truth jacobian:" << std::endl << jacTruth;
  std::cout << "Composite jacobian:" << std::endl << jacComposite2;
  if( !testJacobian( jacComposite2, jacTruth ) )
    {
    std::cout << "Failed getting jacobian for two transforms." << std::endl;
    return EXIT_FAILURE;
    }

//  TODO - test get/set params with only a subset of transforms
//  TODO - test set/clear of transforToOptimzie flags

  return EXIT_SUCCESS;

}