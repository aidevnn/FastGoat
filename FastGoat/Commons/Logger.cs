namespace FastGoat.Commons;

public enum LogLevel
{
    Off,
    Level1,
    Level2
}

public static class Logger
{
    public static LogLevel Level { get; set; } = LogLevel.Off;
}