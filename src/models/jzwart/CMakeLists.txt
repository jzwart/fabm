add_library(fabm_models_jzwart OBJECT
            jzwart_model_library.F90
            doc_decay.F90
           )

add_dependencies(fabm_models_jzwart fabm_base)
