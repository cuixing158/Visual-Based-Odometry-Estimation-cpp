#ifndef _EXTERNAL_H_
#define _EXTERNAL_H_

#ifdef __cplusplus
extern "C" {
#endif

extern void maping_world(int width, int height, unsigned char* image_data, char flag);
extern void init_map_src();

#ifdef __cplusplus
};
#endif

#endif