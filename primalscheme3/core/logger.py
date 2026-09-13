import logging


def close_owned_file_handlers(logger) -> None:
    """Flush and close only file handlers created for this run logger."""

    handlers = tuple(getattr(logger, "_primalscheme_owned_file_handlers", ()))
    for handler in handlers:
        file = None
        try:
            handler.flush()
            file = getattr(getattr(handler, "console", None), "file", None)
            if file is not None:
                file.flush()
        finally:
            logger.removeHandler(handler)
            if file is not None:
                file.close()
            handler.close()
    logger._primalscheme_owned_file_handlers = ()


def setup_rich_logger(logfile: str | None = None):
    from rich.console import Console
    from rich.logging import RichHandler

    console_handler = RichHandler(level=logging.INFO, markup=True, show_path=False)
    logging.basicConfig(
        level="NOTSET",
        format="%(message)s",
        datefmt="[%X]",
        handlers=[console_handler],
    )
    log = logging.getLogger(f"rich-{hash(logfile)}")
    close_owned_file_handlers(log)
    log.setLevel(logging.DEBUG)
    if logfile is not None:
        file_handler = RichHandler(
            level=logging.DEBUG,
            console=Console(file=open(logfile, "w")),
            markup=True,
            show_path=False,
            omit_repeated_times=False,
        )
        file_handler.setFormatter(logging.Formatter("%(message)s", datefmt="[%X]"))
        log.addHandler(file_handler)
        log._primalscheme_owned_file_handlers = (file_handler,)

    return log
