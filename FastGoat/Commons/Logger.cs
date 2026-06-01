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
    public static bool IsOn => Level != LogLevel.Off;
    public static void SetLevel1() => Level = LogLevel.Level1;
    public static void SetLevel2() => Level = LogLevel.Level2;
    public static LogLevel SetOff()
    {
        var lvl = Level;
        Level = LogLevel.Off;
        return lvl;
    }
}