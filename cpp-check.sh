for f in $(find . -name "*.hpp"); do cppcheck --bug-hunting  $f; done
for f in $(find . -name "*.cpp"); do cppcheck --bug-hunting  $f; done
