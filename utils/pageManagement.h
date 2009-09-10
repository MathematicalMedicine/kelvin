#ifndef __PAGE_MANAGEMENT__
#define __PAGE_MANAGEMENT__

void setupSegvHandler (void);
void *allocatePages (int objectSizeInBytes);
void allowReadOnly (void *object, int objectSizeInBytes);
void allowReadWrite (void *object, int objectSizeInBytes);

#endif
