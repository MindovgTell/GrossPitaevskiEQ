add_executable(appcli
    main.cpp
)

target_link_libraries(appcli
    PRIVATE
        GPES          # core library
        # CLI11::CLI11  # CLI parsing
)

# Platform-specific flags for the CLI executable
if(APPLE)
    target_compile_options(appcli
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
    )
elseif(UNIX)
    target_compile_options(appcli
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
    )
endif()