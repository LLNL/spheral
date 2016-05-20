/*--------------------------------------------------------------------------*
 *
 * SpheralDb.c -- Functions that support Spheral's database.
 *
 *--------------------------------------------------------------------------*/

#ifdef USE_POSTGRES

/* Pull in elog() and palloc() functions. You will need postgres.h, 
   util/palloc.h and util/elog.h from the PostGres source tree.  I'm going 
   to bug the PostGreSQL folks about this right now. -JNJ 09.25.2001 */
#include <postgres.h>

/*--------------------------------------------------------------------------*/
/* One dimensional vector. */
typedef struct Vector1d
{
   double x;
} Vector1d;

/* Input function. */
Vector1d*
vector1d_in(char* str)
{
   double x;
   Vector1d* result;
   if (sscanf(str, " (%lf)", &x) != 1) 
   {
      elog(ERROR, "vector1d_in: error in parsing %s", str);
      return NULL;
   } /* end if */
   result = (Vector1d*)palloc(sizeof(Vector1d));
   result->x = x;
   return (result);
}

/* Output function. */
char *
vector1d_out(Vector1d* v)
{
    char *result;
    if (v == NULL)
    {
        return(NULL);
    } /* end if */
    result = (char *) palloc(60);
    sprintf(result, "(%g)", v->x);
    return(result);
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/* Two dimensional vector. */
typedef struct Vector2d
{
   double x, y;
} Vector2d;

/* Input function. */
Vector2d*
vector2d_in(char* str)
{
   double x, y;
   Vector2d* result;
   if (sscanf(str, " (%lf, %lf)", &x, &y) != 2) 
   {
      elog(ERROR, "vector2d_in: error in parsing %s", str);
      return NULL;
   } /* end if */
   result = (Vector2d*)palloc(sizeof(Vector2d));
   result->x = x;
   result->y = y;
   return (result);
}

/* Output function. */
char *
vector2d_out(Vector2d* v)
{
    char *result;
    if (v == NULL)
    {
        return(NULL);
    } /* end if */
    result = (char *) palloc(60);
    sprintf(result, "(%g, %g)", v->x, v->y);
    return(result);
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/* Three dimensional vector. */
typedef struct Vector3d
{
   double x, y, z;
} Vector3d;

/* Input function. */
Vector3d*
vector3d_in(char* str)
{
   double x, y, z;
   Vector3d* result;
   if (sscanf(str, " (%lf, %lf, %lf)", &x, &y, &z) != 3) 
   {
      elog(ERROR, "vector3d_in: error in parsing %s", str);
      return NULL;
   } /* end if */
   result = (Vector3d*)palloc(sizeof(Vector3d));
   result->x = x;
   result->y = y;
   result->z = z;
   return (result);
}

/* Output function. */
char *
vector3d_out(Vector3d* v)
{
    char *result;
    if (v == NULL)
    {
        return(NULL);
    } /* end if */
    result = (char *) palloc(90);
    sprintf(result, "(%g, %g, %g)", v->x, v->y, v->z);
    return(result);
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/* One dimensional tensor. */
typedef struct Tensor1d
{
   double xx;
} Tensor1d;

/* Input function. */
Tensor1d*
tensor1d_in(char* str)
{
   double xx;
   Tensor1d* result;
   if (sscanf(str, " (%lf)", &xx) != 1) 
   {
      elog(ERROR, "tensor1d_in: error in parsing %s", str);
      return NULL;
   } /* end if */
   result = (Tensor1d*)palloc(sizeof(Tensor1d));
   result->xx = xx;
   return (result);
}

/* Output function. */
char *
tensor1d_out(Tensor1d* t)
{
    char *result;
    if (t == NULL)
    {
        return(NULL);
    } /* end if */
    result = (char *) palloc(60);
    sprintf(result, "(%g)", t->xx);
    return(result);
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/* Two dimensional tensor. */
typedef struct Tensor2d
{
   double xx, xy, yx, yy;
} Tensor2d;

/* Input function. */
Tensor2d*
tensor2d_in(char* str)
{
   double xx, xy, yx, yy;
   Tensor2d* result;
   if (sscanf(str, " (%lf, %lf, %lf, %lf)", 
              &xx, &xy, &yx, &yy) != 4) 
   {
      elog(ERROR, "tensor2d_in: error in parsing %s", str);
      return NULL;
   } /* end if */
   result = (Tensor2d*)palloc(sizeof(Tensor2d));
   result->xx = xx;
   result->xy = xy;
   result->yx = yx;
   result->yy = yy;
   return (result);
}

/* Output function. */
char *
tensor2d_out(Tensor2d* t)
{
    char *result;
    if (t == NULL)
    {
        return(NULL);
    } /* end if */
    result = (char *) palloc(120);
    sprintf(result, "(%g, %g, %g, %g)", t->xx, t->xy, t->yx, t->yy);
    return(result);
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/* Three dimensional tensor. */
typedef struct Tensor3d
{
   double xx, xy, xz, yx, yy, yz, zx, zy, zz;
} Tensor3d;

/* Input function. */
Tensor3d*
tensor3d_in(char* str)
{
   double xx, xy, xz, yx, yy, yz, zx, zy, zz;
   Tensor3d* result;
   if (sscanf(str, " (%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf)", 
              &xx, &xy, &xz, &yx, &yy, &yz, &zx, &zy, &zz) != 9) 
   {
      elog(ERROR, "tensor3d_in: error in parsing %s", str);
      return NULL;
   } /* end if */
   result = (Tensor3d*)palloc(sizeof(Tensor3d));
   result->xx = xx;
   result->xy = xy;
   result->xz = xz;
   result->yx = yx;
   result->yy = yy;
   result->yz = yz;
   result->zx = zx;
   result->zy = zy;
   result->zz = zz;
   return (result);
}

/* Output function. */
char *
tensor3d_out(Tensor3d* t)
{
    char *result;
    if (t == NULL)
    {
        return(NULL);
    } /* end if */
    result = (char *) palloc(270);
    sprintf(result, "(%g, %g, %g, %g, %g, %g, %g, %g, %g)", 
            t->xx, t->xy, t->xz, t->yx, t->yy, t->yz, t->zx, t->zy, t->zz);
    return(result);
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/* One dimensional symmetric tensor. */
typedef struct SymTensor1d
{
   double xx;
} SymTensor1d;

/* Input function. */
SymTensor1d*
symtensor1d_in(char* str)
{
   double xx;
   SymTensor1d* result;
   if (sscanf(str, " (%lf)", &xx) != 1) 
   {
      elog(ERROR, "symtensor1d_in: error in parsing %s", str);
      return NULL;
   } /* end if */
   result = (SymTensor1d*)palloc(sizeof(SymTensor1d));
   result->xx = xx;
   return (result);
}

/* Output function. */
char *
symtensor1d_out(SymTensor1d* t)
{
    char *result;
    if (t == NULL)
    {
        return(NULL);
    } /* end if */
    result = (char *) palloc(60);
    sprintf(result, "(%g)", t->xx);
    return(result);
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/* Two dimensional symmetric tensor. */
typedef struct SymTensor2d
{
   double xx, xy, yy;
} SymTensor2d;

/* Input function. */
SymTensor2d*
symtensor2d_in(char* str)
{
   double xx, xy, yy;
   SymTensor2d* result;
   if (sscanf(str, " (%lf, %lf, %lf)", &xx, &xy, &yy) != 3) 
   {
      elog(ERROR, "SymTensor2d_in: error in parsing %s", str);
      return NULL;
   } /* end if */
   result = (SymTensor2d*)palloc(sizeof(SymTensor2d));
   result->xx = xx;
   result->xy = xy;
   result->yy = yy;
   return (result);
}

/* Output function. */
char *
symtensor2d_out(SymTensor2d* t)
{
    char *result;
    if (t == NULL)
    {
        return(NULL);
    } /* end if */
    result = (char *) palloc(90);
    sprintf(result, "(%g, %g, %g)", t->xx, t->xy, t->yy);
    return(result);
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/* Three dimensional symmetric tensor. */
typedef struct SymTensor3d
{
   double xx, xy, xz, yy, yz, zz;
} SymTensor3d;

/* Input function. */
SymTensor3d*
symtensor3d_in(char* str)
{
   double xx, xy, xz, yy, yz, zz;
   SymTensor3d* result;
   if (sscanf(str, " (%lf, %lf, %lf, %lf, %lf, %lf)", 
              &xx, &xy, &xz, &yy, &yz, &zz) != 6) 
   {
      elog(ERROR, "SymTensor3d_in: error in parsing %s", str);
      return NULL;
   } /* end if */
   result = (SymTensor3d*)palloc(sizeof(SymTensor3d));
   result->xx = xx;
   result->xy = xy;
   result->xz = xz;
   result->yy = yy;
   result->yz = yz;
   result->zz = zz;
   return (result);
}

/* Output function. */
char *
symtensor3d_out(SymTensor3d* t)
{
    char *result;
    if (t == NULL)
    {
        return(NULL);
    } /* end if */
    result = (char *) palloc(180);
    sprintf(result, "(%g, %g, %g, %g, %g, %g)", 
                     t->xx, t->xy, t->xz, t->yy, t->yz, t->zz);
    return(result);
}
/*--------------------------------------------------------------------------*/

#endif
