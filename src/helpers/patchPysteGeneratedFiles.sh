for f in $*; do
  echo "Patching " $f
  sed -e 's/std::_Vector_alloc_base/std::vector/g' $f > $f\__
  mv -f $f\__ $f
done
