#~/bin/bash
echo "Running unit tests:"
for i in tests/*_tests
    ./$i_tests
done
echo ""