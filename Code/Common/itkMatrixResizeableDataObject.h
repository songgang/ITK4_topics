/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkMatrixResizeableDataObject_h
#define __itkMatrixResizeableDataObject_h

#include "itkDataObject.h"
#include "itkObjectFactory.h"
#include "vnl/vnl_matrix.h"

namespace itk
{
/**
 * \class MatrixResizeableDataObject
 * \brief Allows for a vnl matrix to be a data
 * object with the flexibility of being resizable.
 *
 * \ingroup DataProcessing
 */

template< typename TItemType >
class MatrixResizeableDataObject:public DataObject, public vnl_matrix< TItemType >
{
public:

  /** Standard class typedefs. */
  typedef MatrixResizeableDataObject Self;
  typedef DataObject                 Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MatrixResizeableDataObject, DataObject);
protected:

  /** Default Constructor. */
  MatrixResizeableDataObject();

  /** Default Destructor. */
  ~MatrixResizeableDataObject();
};
} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_MatrixResizeableDataObject(_, EXPORT, TypeX, TypeY)     \
  namespace itk                                                              \
  {                                                                          \
  _( 1 ( class EXPORT MatrixResizeableDataObject< ITK_TEMPLATE_1 TypeX > ) ) \
  namespace Templates                                                        \
  {                                                                          \
  typedef MatrixResizeableDataObject< ITK_TEMPLATE_1 TypeX >                 \
  MatrixResizeableDataObject##TypeY;                                       \
  }                                                                          \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkMatrixResizeableDataObject+-.h"
#endif

#if ITK_TEMPLATE_TXX
#include "itkMatrixResizeableDataObject.txx"
#endif

#endif
