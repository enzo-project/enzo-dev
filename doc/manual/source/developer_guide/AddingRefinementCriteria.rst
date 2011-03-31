Adding new refinement criteria
==============================


#. Add any new parameters you might need.  (See :ref:`AddingNewParameters`.)
#. Write your code to flag the cells
#. Call your method.

The first point has been discussed elsewhere.

Writing your code to flag cells
-------------------------------

Your code needs to do a couple things:


#. Be named ``FlagCellsToBeRefinedByXXXXXX``, where ``XXXXXX`` is your
   criterion.
#. Increment ``FlaggingField[i]``
#. Count and return the number of flaggged cells.
#. Return -1 on error.

Your code to do the cell flagging can be a grid method.

A minimal code should look like this:

.. code-block:: c

    int grid::FlagCellsToBeRefinedByDensityOverTwelve(){
    
      int NumberOfFlaggedCells = 0;
      for( int i = 0; i< GridDimension[0]*GridDimension[1]*GridDimension[2]; i++ ){
         if( BaryonField[0][i] > 12 ){
           FlaggingField[i] ++;
           NumberOfFlaggedCells ++;
         }
      }
      return NumberOfFlaggedCells;
    
    }

Call your method
----------------

Edit the file ``Grid_SetFlaggingField.C`` In this routine, there's a
loop over the ``CellFlaggingMethod`` array. In this loop, you'll see
code like this:

.. code-block:: c

          /* ==== METHOD 47: By Density over 12 ==== */
    
        case 47:
    
          NumberOfFlaggedCells = this->FlagCellsToBeRefinedByDensityOverTwelve();
          if (NumberOfFlaggedCells < 0) {
            fprintf(stderr, "Error in grid->FlagCellsToBeRefinedByDensityOverTwelve.\n");
            return FAIL;
          }
          break;

So we need to add a few things.


-  Add a new case statement to the switch construct.
-  Set ``NumberOfFlaggedCells`` via the method described above.
-  Don't forget the break; statement.
-  Check for errors.




