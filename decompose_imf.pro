TEMPLATE = subdirs
CONFIG += ordered
SUBDIRS = \
	cpp_utils \ 
	qt_utils \
	decompose_imf_lib \
	decompose_imf_gui \
	convert_matrix \
        decompose_imf_batch \

qt_utils.depends            = cpp_utils
decompose_imf_lib.depends   = cpp_utils qt_utils
decompose_imf_gui.depends   = cpp_utils qt_utils decompose_imf_lib
decompose_imf_batch.depends = cpp_utils qt_utils decompose_imf_lib
convert_matrix.depends      = cpp_utils qt_utils
