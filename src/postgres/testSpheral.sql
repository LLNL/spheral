-----------------------------------------------------------------------------
--
-- testSpheral.sql
--
-- Jeffrey Johnson
-- September 9, 2001
-- Last modified September 26, 2001.
-----------------------------------------------------------------------------
--
-- Test script for spheral++'s database.  Make sure you have done the 
-- following 'fore you try to run this:
--
-- 1. Create the Spheral Database using the createdb command.
-- 2. Create a database user named $USER and endow it with the 
--    priviledge to create databases (make it a "super user").
-- 3. Run the createSpheral.sql script, which will create all the tables 
--    and supporting structures in the database.
--
-- All values will go into the 'test' simulation, and should be deleted 
-- after the tests have finished and displayed.  
-----------------------------------------------------------------------------

-----------------------------------------------------------------------------
--
-- Single value tests
--
-----------------------------------------------------------------------------
-- Write single test values to the value tables.
-- Notice how we have to put quotes around vector and tensor types.
-- Notice also that you need a space between commas separating vector/tensor
-- elements, because I wrote the parsers really lazy-like.  Developers 
-- shouldn't have to deal much with this...
SELECT write_integer('test', 'A', 1);
SELECT write_boolean('test', 'B', true);
SELECT write_string('test', 'C', 'The quick brown fox jumps over lazy dogs');
SELECT write_scalar('test', 'D', 3.14159265);
SELECT write_vector1d('test', 'E', '(1)');
SELECT write_vector2d('test', 'F', '(1, 0)');
SELECT write_vector3d('test', 'G', '(1, 0, 0)');
SELECT write_tensor1d('test', 'H', '(1)');
SELECT write_tensor2d('test', 'I', '(1, 0, 0, 1)');
SELECT write_tensor3d('test', 'J', '(1, 0, 0, 0, 1, 0, 0, 0, 1)');
SELECT write_symtensor1d('test', 'K', '(1)');
SELECT write_symtensor2d('test', 'L', '(1, 0, 1)');
SELECT write_symtensor3d('test', 'M', '(1, 0, 0, 1, 0, 1)');

-- Read the values we've written back out.
SELECT read_integer('test', 'A');
SELECT read_boolean('test', 'B');
SELECT read_string('test', 'C');
SELECT read_scalar('test', 'D');
SELECT read_vector1d('test', 'E');
SELECT read_vector2d('test', 'F');
SELECT read_vector3d('test', 'G');
SELECT read_tensor1d('test', 'H');
SELECT read_tensor2d('test', 'I');
SELECT read_tensor3d('test', 'J');
SELECT read_symtensor1d('test', 'K');
SELECT read_symtensor2d('test', 'L');
SELECT read_symtensor3d('test', 'M');

-- Now delete everything we've added.. :-P
SELECT delete_integer('test', 'A');
SELECT delete_boolean('test', 'B');
SELECT delete_string('test', 'C');
SELECT delete_scalar('test', 'D');
SELECT delete_vector1d('test', 'E');
SELECT delete_vector2d('test', 'F');
SELECT delete_vector3d('test', 'G');
SELECT delete_tensor1d('test', 'H');
SELECT delete_tensor2d('test', 'I');
SELECT delete_tensor3d('test', 'J');
SELECT delete_symtensor1d('test', 'K');
SELECT delete_symtensor2d('test', 'L');
SELECT delete_symtensor3d('test', 'M');

-----------------------------------------------------------------------------
--
-- Field tests
--
-----------------------------------------------------------------------------
-- Write test fields to the field tables.
-- Because of the limitations of PL/PGSQL, we must do much 
-- of this work manually. 
-- Let's use fields that have 10 elements.  This, of course, cannot be 
-- mistaken for an "exhaustive" field test, but we can always augment this 
-- later.  I don't wanna type that much.
-- Notice how we also have to bound field values with quotes.  This means 
-- that fields of vectors/tensors need different quotes around the fields 
-- and the elements.  Ugh...
INSERT INTO scalarfields (tag_id, name_id, elements)
     VALUES (get_tag_id('test'), get_name_id('A'), 
             '{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 
               9.0, 0.0}');
INSERT INTO vector1dfields (tag_id, name_id, elements)
     VALUES (get_tag_id('test'), get_name_id('B'), 
             '{"(1)", "(2)", "(3)", "(4)", "(5)", "(6)", 
               "(7)", "(8)", "(9)", "(0)"}');
INSERT INTO vector2dfields (tag_id, name_id, elements)
     VALUES (get_tag_id('test'), get_name_id('C'), 
             '{"(1, 2)", "(2, 3)", "(3, 4)", "(4, 5)", 
               "(5, 6)", "(6, 7)", "(7, 8)", "(8, 9)", 
               "(9, 0)", "(0, 1)"}');
INSERT INTO vector3dfields (tag_id, name_id, elements)
     VALUES (get_tag_id('test'), get_name_id('D'), 
             '{"(1, 2, 3)", "(2, 3, 4)", "(3, 4, 5)", 
               "(4, 5, 6)", "(5, 6, 7)", "(6, 7, 8)", 
               "(7, 8, 9)", "(8, 9, 0)", "(9, 0, 1)", 
               "(0, 1, 2)"}');
INSERT INTO tensor1dfields (tag_id, name_id, elements)
     VALUES (get_tag_id('test'), get_name_id('E'), 
             '{"(1)", "(2)", "(3)", "(4)", "(5)", "(6)", 
               "(7)", "(8)", "(9)", "(0)"}');
INSERT INTO tensor2dfields (tag_id, name_id, elements)
     VALUES (get_tag_id('test'), get_name_id('F'), 
             '{"(1, 2, 3, 4)", "(2, 3, 4, 5)", "(3, 4, 5, 6)", 
               "(4, 5, 6, 7)", "(5, 6, 7, 8)", "(6, 7, 8, 9)", 
               "(7, 8, 9, 0)", "(8, 9, 0, 1)", "(9, 0, 1, 2)", 
               "(0, 1, 2, 3)"}');
INSERT INTO tensor3dfields (tag_id, name_id, elements)
     VALUES (get_tag_id('test'), get_name_id('G'), 
             '{"(1, 2, 3, 4, 5, 6, 7, 8, 9)", 
               "(2, 3, 4, 5, 6, 7, 8, 9, 0)", 
               "(3, 4, 5, 6, 7, 8, 9, 0, 1)", 
               "(4, 5, 6, 7, 8, 9, 0, 1, 2)", 
               "(5, 6, 7, 8, 9, 0, 1, 2, 3)", 
               "(6, 7, 8, 9, 0, 1, 2, 3, 4)", 
               "(7, 8, 9, 0, 1, 2, 3, 4, 5)", 
               "(8, 9, 0, 1, 2, 3, 4, 5, 6)", 
               "(9, 0, 1, 2, 3, 4, 5, 6, 7)", 
               "(0, 1, 2, 3, 4, 5, 6, 7, 8)"}');
INSERT INTO symtensor1dfields (tag_id, name_id, elements)
     VALUES (get_tag_id('test'), get_name_id('H'), 
             '{"(1)", "(2)", "(3)", "(4)", "(5)", "(6)", 
               "(7)", "(8)", "(9)", "(0)"}');
INSERT INTO symtensor2dfields (tag_id, name_id, elements)
     VALUES (get_tag_id('test'), get_name_id('I'), 
             '{"(1, 2, 3)", "(2, 3, 4)", "(3, 4, 5)", 
               "(4, 5, 6)", "(5, 6, 7)", "(6, 7, 8)", 
               "(7, 8, 9)", "(8, 9, 0)", "(9, 0, 1)", 
               "(0, 1, 2)"}');
INSERT INTO symtensor3dfields (tag_id, name_id, elements)
     VALUES (get_tag_id('test'), get_name_id('J'), 
             '{"(1, 2, 3, 4, 5, 6, 7)", 
               "(2, 3, 4, 5, 6, 7, 8)", 
               "(3, 4, 5, 6, 7, 8, 9)", 
               "(4, 5, 6, 7, 8, 9, 0)", 
               "(5, 6, 7, 8, 9, 0, 1)", 
               "(6, 7, 8, 9, 0, 1, 2)", 
               "(7, 8, 9, 0, 1, 2, 3)", 
               "(8, 9, 0, 1, 2, 3, 4)", 
               "(9, 0, 1, 2, 3, 4, 5)", 
               "(0, 1, 2, 3, 4, 5, 6)"}');

-- Read the field values we've written back out.
-- We also have to do this manually, because of 
-- difficulties in PL/PGSQL of handling arrays.
-- Harumph.
SELECT elements
  FROM scalarfields
 WHERE tag_id = get_tag_id('test')
   AND name_id = get_name_id('A');

SELECT elements
  FROM vector1dfields
 WHERE tag_id = get_tag_id('test')
   AND name_id = get_name_id('B');

SELECT elements
  FROM vector2dfields
 WHERE tag_id = get_tag_id('test')
   AND name_id = get_name_id('C');

SELECT elements
  FROM vector3dfields
 WHERE tag_id = get_tag_id('test')
   AND name_id = get_name_id('D');

SELECT elements
  FROM tensor1dfields
 WHERE tag_id = get_tag_id('test')
   AND name_id = get_name_id('E');

SELECT elements
  FROM tensor2dfields
 WHERE tag_id = get_tag_id('test')
   AND name_id = get_name_id('F');

SELECT elements
  FROM tensor3dfields
 WHERE tag_id = get_tag_id('test')
   AND name_id = get_name_id('G');

SELECT elements
  FROM symtensor1dfields
 WHERE tag_id = get_tag_id('test')
   AND name_id = get_name_id('H');

SELECT elements
  FROM symtensor2dfields
 WHERE tag_id = get_tag_id('test')
   AND name_id = get_name_id('I');

SELECT elements
  FROM symtensor3dfields
 WHERE tag_id = get_tag_id('test')
   AND name_id = get_name_id('J');

-- Now delete everything we've added.. :-P
SELECT delete_scalarfield('test', 'A');
SELECT delete_vector1dfield('test', 'B');
SELECT delete_vector2dfield('test', 'C');
SELECT delete_vector3dfield('test', 'D');
SELECT delete_tensor1dfield('test', 'E');
SELECT delete_tensor2dfield('test', 'F');
SELECT delete_tensor3dfield('test', 'G');
SELECT delete_symtensor1dfield('test', 'H');
SELECT delete_symtensor2dfield('test', 'I');
SELECT delete_symtensor3dfield('test', 'J');

