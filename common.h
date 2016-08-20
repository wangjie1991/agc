#ifndef _AUDIO_GOLDENC_COMMON_H_
#define _AUDIO_GOLDENC_COMMON_H_

template <typename Type>
static inline void Malloc(Type** data, size_t num_elems) {
    *data = new Type[num_elems];
    if (*data == NULL) {
        std::cerr << "OOM" << std::endl;
        exit(1);
    }
}

#endif
